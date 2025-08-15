#that was for me to make the R compendium
#let's make my R compendium
#remotes::install_deps()
#https://frbcesab.github.io/rcompendium/articles/rcompendium.html
#I sometime got issue loading new pacakge from GitHub in case
#https://community.rstudio.com/t/unable-to-install-packages-from-github/124372
#Sys.unsetenv("GITHUB_PAT")

#library(rcompendium)
#library(mgcv)
#add_description(given = 'Valentin', family = 'Journé',
#                email = 'journe.valentin@gmail.com')
#add_license('CC BY 4.0')
#add_readme_rmd(given = 'Valentin', family = 'Journé')
# add_code_of_conduct(email = 'journe.valentin@gmail.com')
# add_citation(
#   given = 'Valentin',
#   family = 'Journé',
#   organisation = 'Department of Biology, Faculty of Science, Kyushu University, Fukuoka, Japan',
#   open = TRUE,
#   overwrite = FALSE,
#   quiet = FALSE
# )
# add_dependencies('.')

#PACKAGE MANAGEMENT
#to dev tool
# devtools::load_all(here::here())
# devtools::document()
# devtools::check()
# devtools::install()
# report <- pkgnet::CreatePackageReport(
#   pkg_name = "weatheRcues"
# )
# browseURL(report)

#logo
# sysfonts::font_files() %>% View()
# sysfonts::font_add("Futura", "Futura.ttc")
# p = magick::image_read(
#   "/Users/valentinjourne/Documents/REDACTION/MASTREE_detection_cues/logo.png"
# )
# hexSticker::sticker(
#   subplot = p,
#   package = "weatheRcues",
#   p_size = 14,
#   p_y = 1.6,
#   p_color = "#FDFEFE", #"#FDF6E3",
#   s_x = 1,
#   s_y = 0.75,
#   s_width = 1,
#   s_height = 1,
#   h_color = "#FDCB6E", #"#FFC857",
#   h_fill = "#5E4FA2", #"#CA6702",
#   h_size = 1,
#   p_family = "Futura",
#   filename = "man/figures/logo2.png"
# )
