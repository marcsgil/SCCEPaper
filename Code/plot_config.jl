using CairoMakie, MathTeXEngine, ColorSchemes

colors = colorschemes[:seaborn_bright]
theme = Theme(
    fontsize=32,
    fonts = (; regular = texfont()),
    Axis=(xlabelsize=36, xlabelpadding=0,
    ylabelsize=36, ylabelpadding=0,
    xticklabelsize = 28, yticklabelsize=28,
    xgridvisible=false, ygridvisible=false,
    xtickalign=1, ytickalign=1,
    yticksize=12, xticksize=12),
    Lines = (linewidth = 5,),
    Scatter = (markersize=20,),
    Legend = (labelsize=28,)
)

set_theme!()
set_theme!(theme)