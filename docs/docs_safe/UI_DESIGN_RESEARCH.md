# UI Design Research for DeconVoVo

Research date: 2026-03-25

---

## Table of Contents

1. [PySide6/Qt6 Window Transparency on Windows](#1-pyside6qt6-window-transparency-on-windows)
2. [Beautiful Scientific Desktop App UI Examples](#2-beautiful-scientific-desktop-app-ui-examples)
3. [UI Design Principles for Scientific Tools](#3-ui-design-principles-for-scientific-tools)
4. [Actionable Recommendations for DeconVoVo](#4-actionable-recommendations-for-deconvovo)

---

## 1. PySide6/Qt6 Window Transparency on Windows

### 1.1 Why Our Current Approach Only Partially Works

Our current code in `gui/app.py` uses three methods:

```python
self.setAttribute(Qt.WA_TranslucentBackground, True)
dwm.DwmEnableBlurBehindWindow(hwnd, ...)
dwm.DwmExtendFrameIntoClientArea(hwnd, MARGINS(-1,-1,-1,-1))
dwm.DwmSetWindowAttribute(hwnd, 38, c_int(3), 4)  # DWMWA_SYSTEMBACKDROP_TYPE = Acrylic
```

**The problem:** `DwmEnableBlurBehindWindow` is a legacy Win7 Aero API. On Windows
10/11, it only blurs the DWM-rendered frame area (the titlebar and edges), NOT the
client area. `DwmExtendFrameIntoClientArea` with MARGINS(-1,-1,-1,-1) extends the
glass frame into the client area, but on Win10+ this just makes the extended area
black/transparent -- it does NOT blur.

The `DWMWA_SYSTEMBACKDROP_TYPE = 3` (Acrylic) call is the correct Win11 22H2+ API,
but it ALSO only applies the backdrop behind the non-client area (titlebar + frame).
Per Microsoft docs: "Draw the backdrop material effect corresponding to a transient
window behind the entire window bounds" -- but in practice for Win32 apps, Windows
only renders the backdrop behind the non-client region unless the window explicitly
opts in via WinUI/XAML composition.

### 1.2 The Three Methods for Window Blur/Transparency

#### Method A: SetWindowCompositionAttribute (Undocumented API)

This is the only method that provides FULL CLIENT AREA blur/acrylic on Win10/11
in a Win32 (non-UWP) application. It is used by virtually all Python blur libraries.

**How it works:**

```python
import ctypes
from ctypes.wintypes import DWORD, BOOL, HWND

user32 = ctypes.windll.user32
dwm = ctypes.windll.dwmapi

class ACCENTPOLICY(ctypes.Structure):
    _fields_ = [
        ("AccentState",   ctypes.c_uint),  # Key field
        ("AccentFlags",   ctypes.c_uint),
        ("GradientColor", ctypes.c_uint),  # ABGR format
        ("AnimationId",   ctypes.c_uint),
    ]

class WINDOWCOMPOSITIONATTRIBDATA(ctypes.Structure):
    _fields_ = [
        ("Attribute",  ctypes.c_int),
        ("Data",       ctypes.POINTER(ctypes.c_int)),
        ("SizeOfData", ctypes.c_size_t),
    ]

# AccentState values:
# 0 = ACCENT_DISABLED
# 1 = ACCENT_ENABLE_GRADIENT
# 2 = ACCENT_ENABLE_TRANSPARENTGRADIENT
# 3 = ACCENT_ENABLE_BLURBEHIND         (Win10 standard blur)
# 4 = ACCENT_ENABLE_ACRYLICBLURBEHIND  (Win10 1803+ acrylic)
# 5 = ACCENT_ENABLE_HOSTBACKDROP       (Win11 Mica-like)

WCA_ACCENT_POLICY = 19
WCA_USEDARKMODECOLORS = 26

def enable_acrylic(hwnd, gradient_color=0x40121212, dark=True):
    """Apply full-window acrylic blur. gradient_color is ABGR hex."""
    accent = ACCENTPOLICY()
    accent.AccentState = 4  # ACCENT_ENABLE_ACRYLICBLURBEHIND
    accent.AccentFlags = 2
    accent.GradientColor = gradient_color

    data = WINDOWCOMPOSITIONATTRIBDATA()
    data.Attribute = WCA_ACCENT_POLICY
    data.SizeOfData = ctypes.sizeof(accent)
    data.Data = ctypes.cast(ctypes.pointer(accent), ctypes.POINTER(ctypes.c_int))

    user32.SetWindowCompositionAttribute(int(hwnd), data)

    if dark:
        data.Attribute = WCA_USEDARKMODECOLORS
        user32.SetWindowCompositionAttribute(int(hwnd), data)
```

**Requirements on the Qt side:**
```python
self.setAttribute(Qt.WA_TranslucentBackground, True)
self.setStyleSheet("background-color: rgba(0, 0, 0, 0)")
```

**Pros:**
- Works on Windows 10 1803+ and Windows 11
- Applies blur/acrylic to the ENTIRE window including client area
- Tint color is controllable via GradientColor (ABGR format)
- No need for frameless window

**Cons:**
- **Undocumented API** -- could break in future Windows updates
- **Window dragging is extremely laggy** on Win10 (known bug, no fix)
- Win10 1903+ introduced additional lag when resizing
- AccentState=5 (HostBackdrop/Mica-like) is undocumented and unreliable

#### Method B: DwmSetWindowAttribute with DWMWA_SYSTEMBACKDROP_TYPE (Official Win11 API)

This is the **official, documented** Microsoft API but ONLY works on Windows 11 22H2+
(Build 22621).

```python
DWMWA_SYSTEMBACKDROP_TYPE = 38

# Backdrop types:
# DWMSBT_AUTO            = 0  (let DWM decide)
# DWMSBT_NONE            = 1  (no backdrop)
# DWMSBT_MAINWINDOW      = 2  (Mica)
# DWMSBT_TRANSIENTWINDOW = 3  (Acrylic)
# DWMSBT_TABBEDWINDOW    = 4  (Mica Alt)

def enable_mica(hwnd, dark=True):
    dwm = ctypes.windll.dwmapi
    # Enable dark mode first
    DWMWA_USE_IMMERSIVE_DARK_MODE = 20
    val = ctypes.c_int(1 if dark else 0)
    dwm.DwmSetWindowAttribute(hwnd, DWMWA_USE_IMMERSIVE_DARK_MODE,
                               ctypes.byref(val), ctypes.sizeof(val))
    # Apply Mica backdrop
    backdrop = ctypes.c_int(2)  # DWMSBT_MAINWINDOW = Mica
    dwm.DwmSetWindowAttribute(hwnd, DWMWA_SYSTEMBACKDROP_TYPE,
                               ctypes.byref(backdrop), ctypes.sizeof(backdrop))
```

**Critical limitation:** For Win32 apps, this backdrop is rendered **behind the
non-client area only** (titlebar). Getting it behind the client area requires either:
1. Using a frameless window and extending the frame, OR
2. Using WinUI 3 / XAML Islands (not available in PySide6)

The `PyQt-Frameless-Window` library by zhiyiYo solves this by removing the standard
window frame and reimplementing the titlebar, then applying the Mica backdrop to the
entire frameless surface.

**Pros:**
- Official, documented, stable API
- Best visual quality (Mica adapts to wallpaper)
- No lag issues

**Cons:**
- Windows 11 22H2+ only
- Requires frameless window approach for full-window effect
- More complex implementation

#### Method C: DwmEnableBlurBehindWindow (Legacy, DO NOT USE)

This is the Win7 Aero Glass API. On Windows 10/11 it produces a faint, ugly blur
effect only on the extended frame region. **Not recommended for modern apps.**

### 1.3 GitHub Projects for Window Transparency

| Project | Stars | URL | Method | PySide6 | Notes |
|---------|-------|-----|--------|---------|-------|
| **PyQt-Fluent-Widgets** (zhiyiYo) | 7,700 | [GitHub](https://github.com/zhiyiYo/PyQt-Fluent-Widgets) | Acrylic + Mica via PyQt-Frameless-Window | Yes (PySide6 branch) | Best-in-class Fluent Design. GPLv3 for non-commercial. Has AcrylicLabel, navigation, settings cards. |
| **PyQt-Frameless-Window** (zhiyiYo) | 730 | [GitHub](https://github.com/zhiyiYo/PyQt-Frameless-Window) | Win10 Acrylic (SWCA), Win11 Mica (DwmSetWindowAttribute), Win7 Aero, macOS blur | Yes (PySide6 branch) | The foundation for PyQt-Fluent-Widgets. Frameless + custom titlebar + blur. |
| **PythonBlurBehind** (Peticali) | 66 | [GitHub](https://github.com/Peticali/PythonBlurBehind) | SetWindowCompositionAttribute with ACCENT_ENABLE_ACRYLICBLURBEHIND | PySide2 examples only (works with any HWND) | Standalone blur. Simple API: `GlobalBlur(hwnd, Acrylic=True, Dark=True)`. Code needs modernization per author. |
| **py-window-styles** (Akascape) | 576 | [GitHub](https://github.com/Akascape/py-window-styles) | Multiple: Mica, Acrylic, Aero, Transparent | Any framework (Tkinter, PyQt, PySide, etc.) | Broadest framework support. `pywinstyles.apply_style(window, "mica")`. CC0 license. |
| **BlurWindow** (PyPI) | ~50 | [PyPI](https://pypi.org/project/BlurWindow/) | SetWindowCompositionAttribute | PySide/Tkinter | Similar to PythonBlurBehind. |
| **CustomQt** (ultrasploit) | 3 | [GitHub](https://github.com/ultrasploit/CustomQt) | DWM shadows/accents/blur, WM_NCHITTEST | PySide6 native | Lightweight. Custom titlebars + optional acrylic. Very new/small. |
| **winmica** (amnweb) | 4 | [GitHub](https://github.com/amnweb/winmica) | DwmSetWindowAttribute DWMWA_SYSTEMBACKDROP_TYPE | PyQt6 (should work with PySide6) | Simple `EnableMica(hwnd, BackdropType.MICA)` API. Win11 22H2+ only. |
| **framelesshelper** (wangwenx190) | ~2,000 | [GitHub](https://github.com/wangwenx190/framelesshelper) | Native DWM blur for Win7-11, macOS, Linux | C++ Qt (Python via bindings) | Moved to qwindowkit. Most technically complete. |

### 1.4 The Correct Approach for DeconVoVo

**Recommendation: Use `SetWindowCompositionAttribute` with `ACCENT_ENABLE_BLURBEHIND`
(AccentState=3) for a subtle blur, NOT full acrylic (AccentState=4).**

Why:
- AccentState=3 (standard blur) does NOT have the window-dragging lag issue
- AccentState=4 (acrylic) causes severe lag on Win10 and some Win11 configs
- The blur still looks modern and professional with a semi-transparent tinted overlay
- Combined with `WA_TranslucentBackground` and rgba background colors in QSS, this
  gives a frosted-glass appearance

**Alternative: Use `py-window-styles` as a thin wrapper:**
```python
import pywinstyles
pywinstyles.apply_style(window, "acrylic")  # or "mica" on Win11
```

**Best alternative: If we want rock-solid reliability, skip real OS-level
transparency entirely.** Use a semi-transparent dark background with carefully chosen
colors (e.g., `rgba(26, 26, 36, 240)`) and subtle gradient overlays. Many
"modern-looking" apps (PyDracula, PyOneDark) use NO real transparency and still look
excellent. The visual difference is minimal when the app is in the foreground.

---

## 2. Beautiful Scientific Desktop App UI Examples

### 2.1 Modern UI Framework Projects

#### PyQt-Fluent-Widgets (zhiyiYo) -- 7,700 stars
- **URL:** https://github.com/zhiyiYo/PyQt-Fluent-Widgets
- **Look:** Microsoft Fluent Design System (WinUI 3 style)
- **Key features:** Navigation sidebar with EXPAND/COMPACT/MENU modes, SettingCards
  (switch, range, color, combobox, hyperlink), AcrylicLabel, InfoBar notifications,
  smooth scrolling, Pivot and TabBar navigation, dialog system
- **Sidebar:** `NavigationInterface` with `NavigationPanel` -- three zones (top,
  scroll, bottom). Icons + text that collapse to icon-only. Responsive breakpoints.
- **Themes:** Dark/light with Mica/Acrylic backdrop integration
- **Visual style:** Clean, minimal, generous padding. Looks like a native Windows 11 app.
- **License:** GPLv3 (free for non-commercial)

#### Modern_GUI_PyDracula (Wanderson-Magalhaes) -- 3,000 stars
- **URL:** https://github.com/Wanderson-Magalhaes/Modern_GUI_PyDracula_PySide6_or_PyQt6
- **Look:** Dark theme with Dracula color palette (purple/pink accents on dark bg)
- **Key features:** Sidebar navigation (icon + text, collapsible), custom titlebar,
  theme switching (dark/light), modular page system, QSS-based theming
- **Design pattern:** Frameless window with custom minimize/maximize/close. Left sidebar
  with navigation icons. Content area with stacked pages. Top bar with breadcrumbs.
- **Architecture:** `main.py` + `modules/` (logic) + `widgets/` (custom) + `themes/` (QSS)
- **License:** MIT

#### PyOneDark (Wanderson-Magalhaes) -- 1,100 stars
- **URL:** https://github.com/Wanderson-Magalhaes/PyOneDark_Qt_Widgets_Modern_GUI
- **Look:** Atom One Dark theme inspired. Darker and more muted than PyDracula.
- **Key features:** Custom left menu with icons, tooltip-based navigation, custom
  widgets system, settings via JSON, Qt Designer integration
- **Design pattern:** Similar to PyDracula but with more widget customization. Two-tier
  sidebar (thin icon bar + expandable panel).
- **License:** MIT with donation request for commercial use

#### qt-material (UN-GCPDS) -- 2,800 stars
- **URL:** https://github.com/UN-GCPDS/qt-material
- **Look:** Google Material Design adaptation. Bold accent colors on dark/light.
- **Key features:** 19 built-in themes (9 dark, 10 light), custom XML theme definition,
  runtime theme switching, density scaling (-2 to +2), theme export to standalone QSS,
  accent buttons (danger/warning/success), interactive theme editor
- **Strengths for scientific apps:** Density scaling lets you pack more controls.
  Accent status colors (danger=red, warning=orange, success=green) useful for
  validation feedback.
- **License:** BSD-2-Clause

#### QDarkStyleSheet (ColinDuquesnoy) -- 3,100 stars
- **URL:** https://github.com/ColinDuquesnoy/QDarkStyleSheet
- **Look:** Balanced dark/light themes developed with the Spyder IDE team
- **Key features:** SCSS-based palette system, comprehensive widget coverage,
  designed for scientific IDE (Spyder), PyQtGraph compatible
- **Why it matters:** Designed WITH a scientific tool (Spyder). Proven in production.
- **License:** MIT (code), CC-BY (images)

#### PyQtDarkTheme (5yutan5) -- 741 stars
- **URL:** https://github.com/5yutan5/PyQtDarkTheme
- **Look:** Flat, modern dark/light. Syncs with OS theme and accent colors.
- **Key features:** Auto OS theme sync (Win/Mac/Linux), flat design, fixes Qt version
  styling inconsistencies, custom QPalette, sidebar property for QToolbar
- **License:** MIT

#### QtModernRedux (robertkist) -- 36 stars
- **URL:** https://github.com/robertkist/qtmodernredux
- **Look:** Modern frameless with custom titlebar and drop shadows
- **Key features:** Custom window chrome, DPI-aware, title-area widgets (like Chrome
  tabs), improved QTableView/QListView styling, cross-platform consistency
- **License:** MIT

#### 24-Modern-Desktop-GUI (KhamisiKibet) -- 107 stars
- **URL:** https://github.com/KhamisiKibet/24-Modern-Desktop-GUI
- **Look:** Professional with dark/light themes, custom widgets
- **Key features:** Tutorial-driven, Qt Designer workflow, JSON config system,
  theme hot-reload, PyInstaller packaging
- **License:** Not specified

### 2.2 Scientific Visualization Tools

#### PyQtGraph -- 4,300 stars
- **URL:** https://github.com/pyqtgraph/pyqtgraph
- **Purpose:** Fast 2D/3D scientific plotting for PySide6/PyQt6
- **Features:** Real-time data display, image analysis, ROI tools, export to
  matplotlib/HDF5, GPU-accelerated OpenGL 3D
- **Used by:** Orange3, ACQ4, neurophysiology tools
- **Dark theme compatibility:** Works with QDarkStyleSheet and qt-material

#### Matplotlib in PySide6
- **Tutorial:** https://www.pythonguis.com/tutorials/pyside6-plotting-matplotlib/
- **Approach:** Embed `FigureCanvasQTAgg` widget in PySide6 layouts
- **Our current approach:** Already using matplotlib for plots

### 2.3 Design Patterns Observed Across Projects

**Sidebar Navigation (universal pattern):**
- Left sidebar, 48-60px wide when collapsed (icons only), 200-280px expanded
- Three zones: top (brand/main nav), middle (scrollable pages), bottom (settings/user)
- Icons are 20-24px, with 8-12px padding
- Active item highlighted with accent color bar (2-3px left border) or filled background
- Hover state: subtle background highlight

**Content Area:**
- Pages in a QStackedWidget, switched by sidebar selection
- Consistent padding: 16-24px on all sides
- Card-based grouping for related controls
- Headers with clear hierarchy (title > subtitle > body)

**Color Schemes (dark theme):**
- Background: #1a1a24 to #2b2b3d (very dark blue-gray)
- Surface/Card: +10-15% lighter than background
- Border: subtle, 1px, ~20% opacity white
- Text primary: #e0e0e0 to #ffffff
- Text secondary: #a0a0b0 (60-70% opacity)
- Accent: one strong color (teal, blue, purple) used sparingly
- Status: red=error, amber=warning, green=success, blue=info

---

## 3. UI Design Principles for Scientific Tools

### 3.1 Spacing and Layout

**Grid System:**
- Use a consistent base unit (8px is standard). All spacing should be multiples: 8, 16, 24, 32, 40px.
- This creates visual rhythm and subconscious order.

**Padding Guidelines:**
- Window margins: 16-24px
- Card/group internal padding: 12-16px
- Between related controls: 8px
- Between control groups: 16-24px
- Between sections: 32px

**Alignment:**
- Left-align labels and controls consistently
- Right-align numeric data in tables
- Center-align status indicators and icons
- Use consistent label widths (form labels same width via QFormLayout or grid)

### 3.2 Typography Hierarchy

**Font Size Scale (for 96 DPI):**
- Window title: 16-18pt, bold, accent color
- Section title: 13-14pt, semibold
- Subsection/card title: 11-12pt, semibold
- Body text / labels: 10pt, regular
- Secondary/caption text: 9pt, regular, dimmed color
- Table data: 9-10pt, regular (monospace for numbers recommended)

**Key rules:**
- Maximum 2-3 font sizes per view
- Use weight (bold/semibold) and color to create hierarchy, not just size
- Monospace font for numeric data, paths, and scientific notation
- Never use more than 2 font families

### 3.3 Color Usage

**Accent Color Strategy:**
- Pick ONE primary accent (we use teal/cyan: #00bfa5)
- Use it for: active nav items, primary buttons, focused inputs, links, key metrics
- Use desaturated/dimmed variant for secondary emphasis
- Never use accent for large background areas

**Status Colors:**
- Error/failure: #ef5350 (red)
- Warning/caution: #ffa726 (amber)
- Success/complete: #66bb6a (green)
- Info/neutral: #42a5f5 (blue)
- Always pair color with icon or text (accessibility)

**Surface Hierarchy (dark theme):**
- Layer 0 (window bg): darkest
- Layer 1 (sidebar, cards): slightly lighter
- Layer 2 (input fields, nested cards): slightly lighter again
- Layer 3 (hover/active states): lightest

### 3.4 Information Density vs. Whitespace

**The Scientific App Challenge:**
Scientists need to see many parameters, results, and plots simultaneously.
"Dashboard whitespace" wastes their screen real estate. But dense UIs without
structure look chaotic and unprofessional.

**Solution -- Structured Density:**
- Group related controls visually with cards or subtle background regions
- Use horizontal dividers (1px, subtle) between logical sections
- Compact spacing WITHIN groups, generous spacing BETWEEN groups
- Collapsible sections for "advanced" parameters
- Tables with alternating row colors for scanability
- Status bars and info strips at top/bottom for persistent context

**Information Priority:**
1. Primary data visualization (plots, heatmaps) -- largest area, most prominent
2. Current results / metrics -- visible at a glance
3. Input parameters -- accessible but not dominant
4. Secondary controls / settings -- available but tucked away
5. Status / log / metadata -- persistent but minimal

### 3.5 Control Grouping and Visual Flow

**Logical Grouping:**
- Group by workflow step (input > process > output)
- Each group gets a card or bordered section with a clear title
- Related controls within a group should be visually connected

**Visual Flow:**
- Left to right, top to bottom (Western reading order)
- Sidebar (navigate) -> Content (work) -> Results (see)
- Within content: Parameters (top) -> Action button -> Results (bottom)
- Eye should naturally flow from "what do I set?" to "what do I get?"

**Progressive Disclosure:**
- Show essential controls by default
- Hide advanced/rarely-used options behind expandable sections or "Advanced..." buttons
- Use tabs for parallel content that shares the same space
- Tooltips for parameter descriptions (not inline text)

### 3.6 Scientific Software-Specific Principles

From research on usability in scientific software:

1. **Precision over simplicity:** Scientists need exact control. Numeric spinboxes
   with editable text, not just sliders. Show units always.

2. **Reproducibility support:** Show parameter state clearly. Allow saving/loading
   parameter sets. Display file paths and versions.

3. **Feedback at every step:** Progress bars for long operations. Status messages
   for each pipeline stage. Clear error messages with actionable suggestions.

4. **Separation of concerns:** UI code separate from computation code. Parameters
   stored in files/configs, not hidden in GUI state.

5. **Command-line parity:** Anything the GUI does should be reproducible from CLI.
   Show the equivalent command or parameter file.

6. **Domain vocabulary:** Use terms the scientist knows (m/z, drift time, CCS),
   not programmer terms (array index, buffer size).

7. **Metadata visibility:** File names, sample IDs, timestamps, instrument
   parameters should be readily accessible, not buried in menus.

---

## 4. Actionable Recommendations for DeconVoVo

### 4.1 Transparency Decision

**Recommended approach:** Drop real OS-level blur/acrylic transparency. Instead:
- Use solid dark backgrounds with carefully chosen colors to simulate depth
- Use `rgba()` backgrounds in QSS for internal layering (cards on surface)
- Add subtle 1px borders and box-shadows (via QGraphicsDropShadowEffect) for depth
- This avoids all the DWM API headaches, works cross-platform, and looks just as good

If we really want the effect:
- Use `py-window-styles` library (`pywinstyles.apply_style(window, "acrylic")`)
  for Win10/11 with graceful fallback
- Accept the dragging lag on Win10 as a tradeoff
- Set AccentState=3 (standard blur, no lag) instead of 4 (acrylic, laggy)

### 4.2 Theme Library Decision

**Best options ranked:**

1. **qt-material** (2,800 stars, BSD-2) -- Best for us. Drop-in stylesheet, 19 themes,
   density scaling for scientific UIs, status accent colors, no structural changes needed.
   `apply_stylesheet(app, theme='dark_teal.xml')` and done.

2. **PyQtDarkTheme** (741 stars, MIT) -- Simplest. One line: `qdarktheme.setup_theme("dark")`.
   Auto-syncs with OS. Sidebar support via QToolbar property.

3. **QDarkStyleSheet** (3,100 stars, MIT) -- Proven in Spyder. More conservative look.
   SCSS-based palette for customization.

4. **Custom QSS** (what we have now) -- Maximum control, most maintenance. Keep if we
   want pixel-perfect brand identity.

5. **PyQt-Fluent-Widgets** (7,700 stars, GPLv3) -- Most beautiful, but GPLv3 is viral.
   Commercial license available. Overkill unless we need the full navigation system.

### 4.3 Navigation Pattern

Current DeconVoVo uses a left sidebar with icon-text buttons. Improvements:
- Add visual indicator for active page (accent-colored left border, 3px)
- Add hover state (subtle background highlight)
- Add section dividers in sidebar between logical groups
- Consider collapsible sidebar (icons-only mode) for more content area

### 4.4 Specific Improvements to Implement

**Quick wins (no library changes):**
1. Increase card corner radius to 8px (currently likely 4px or 0)
2. Add consistent 8px-grid spacing throughout
3. Add alternating row colors in all QTableViews
4. Use monospace font for numeric data columns
5. Add subtle QGraphicsDropShadowEffect to cards/panels
6. Increase contrast between sidebar and content background
7. Add active-state indicator (accent border) to sidebar items

**Medium effort:**
1. Switch to `qt-material` or `PyQtDarkTheme` for consistent widget styling
2. Add status color coding (green=complete, amber=processing, red=error)
3. Add progress feedback for pipeline operations
4. Group parameters into collapsible cards with clear titles
5. Add icons to sidebar items (use Qt built-in or Material Design icons)

**Larger effort:**
1. Implement responsive sidebar (collapsible to icons on narrow windows)
2. Add parameter save/load (JSON presets)
3. Add theme toggle that persists across sessions (already started)
4. Consider embedded PyQtGraph for real-time interactive plots
5. Add "pipeline status" strip showing current stage and progress

### 4.5 Color Palette Recommendation

For our dark theme with teal accent (#00bfa5):

```
Background (Layer 0):   #1a1a2e   (deep navy-black)
Surface (Layer 1):      #222236   (sidebar, cards)
Input fields:           #2a2a40   (slightly lighter)
Hover:                  #32324e   (noticeable but subtle)
Border:                 #3a3a55   (visible but not harsh)
Text primary:           #e8e8f0   (near white, slight blue tint)
Text secondary:         #8888a0   (readable but receded)
Text disabled:          #555570   (clearly muted)
Accent:                 #00bfa5   (teal -- active items, links, primary buttons)
Accent dim:             #00897b   (secondary emphasis)
Error:                  #ef5350
Warning:                #ffa726
Success:                #66bb6a
Info:                   #42a5f5
```

---

## Appendix: All GitHub Repos Referenced

| Project | Stars | License | URL |
|---------|-------|---------|-----|
| PyQt-Fluent-Widgets | 7,700 | GPLv3 | https://github.com/zhiyiYo/PyQt-Fluent-Widgets |
| PyQtGraph | 4,300 | MIT | https://github.com/pyqtgraph/pyqtgraph |
| QDarkStyleSheet | 3,100 | MIT | https://github.com/ColinDuquesnoy/QDarkStyleSheet |
| PyDracula Modern GUI | 3,000 | MIT | https://github.com/Wanderson-Magalhaes/Modern_GUI_PyDracula_PySide6_or_PyQt6 |
| qt-material | 2,800 | BSD-2 | https://github.com/UN-GCPDS/qt-material |
| framelesshelper | ~2,000 | MIT | https://github.com/wangwenx190/framelesshelper |
| PyOneDark | 1,100 | MIT | https://github.com/Wanderson-Magalhaes/PyOneDark_Qt_Widgets_Modern_GUI |
| PyQtDarkTheme | 741 | MIT | https://github.com/5yutan5/PyQtDarkTheme |
| PyQt-Frameless-Window | 730 | MIT | https://github.com/zhiyiYo/PyQt-Frameless-Window |
| py-window-styles | 576 | CC0 | https://github.com/Akascape/py-window-styles |
| 24-Modern-Desktop-GUI | 107 | -- | https://github.com/KhamisiKibet/24-Modern-Desktop-GUI |
| PythonBlurBehind | 66 | MIT | https://github.com/Peticali/PythonBlurBehind |
| QtModernRedux | 36 | MIT | https://github.com/robertkist/qtmodernredux |
| winmica | 4 | MIT | https://github.com/amnweb/winmica |
| CustomQt | 3 | MIT | https://github.com/ultrasploit/CustomQt |

## Appendix: Key Web Resources

- [DWM_SYSTEMBACKDROP_TYPE (Microsoft Docs)](https://learn.microsoft.com/en-us/windows/win32/api/dwmapi/ne-dwmapi-dwm_systembackdrop_type)
- [DWMWINDOWATTRIBUTE (Microsoft Docs)](https://learn.microsoft.com/en-us/windows/win32/api/dwmapi/ne-dwmapi-dwmwindowattribute)
- [Mica Material Design Guide (Microsoft)](https://learn.microsoft.com/en-us/windows/apps/design/style/mica)
- [Acrylic Material Design Guide (Microsoft)](https://learn.microsoft.com/en-us/windows/apps/design/style/acrylic)
- [PySide6 Matplotlib Embedding Tutorial](https://www.pythonguis.com/tutorials/pyside6-plotting-matplotlib/)
- [PySide6 PyQtGraph Tutorial](https://www.pythonguis.com/tutorials/pyside6-plotting-pyqtgraph/)
- [PySide6 QTableView with Pandas](https://www.pythonguis.com/tutorials/pyside6-qtableview-modelviews-numpy-pandas/)
- [Good Usability Practices in Scientific Software (arXiv)](https://arxiv.org/pdf/1709.00111)
- [Qt-Material Documentation](https://qt-material.readthedocs.io/)
- [PyQt-Fluent-Widgets Navigation Docs](https://pyqt-fluent-widgets.readthedocs.io/en/latest/navigation.html)
- [UI Design Principles (Figma)](https://www.figma.com/resource-library/ui-design-principles/)
- [Typography in UX/UI Guide](https://supercharge.design/blog/typography-in-ux-ui-a-complete-guide)
- [Data Visualization UI Best Practices](https://www.transcenda.com/insights/data-visualization-ui-best-practices-and-winning-approaches)
