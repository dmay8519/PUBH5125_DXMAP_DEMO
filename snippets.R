```{r GeoDa, eval=FALSE}

geoda <- utils::read.csv(text = c("Variable,Description",
                                  "CODENO,Code converted to numeric (drop w prefix)",
                                  "AREA,District polygon area",
                                  "PERIMETER,District polygon perimeter",
                                  "RECORD_ID,Unique ID",
                                  "DISTRICT,District number 1-56",
                                  "NAME,Name of districts from Cressie (1993)",
                                  "CODE,District code from WinBugs",
                                  "CANCER,Lip cancer cases from Cressie (1993)",
                                  "POP,Population years at risk from Cressie (1993)",
                                  "CEXP,Expected cases from Lawson et al. (1999)",
                                  "AFF,Outdoor industry from Lawson et al. (1999)")) |>
  dplyr::mutate(`Corrected Description` = base::gsub("Lawson et al\\.", "Stern & Cressie", Description))

even <- geoda %>%
  dplyr::mutate(lines = 1:dplyr::n()) %>%
  dplyr::filter((lines - 1) %% 2 == 1) %>%
  dplyr::pull(lines)

# 1:base::nrow(geoda)[1:base]
#   
# 
# foo <- knitr::kable(geoda,
#              align = c("l", "l", "l"), 
#              row.names = FALSE,
#              escape = FALSE) %>%
# 
# 
# 
#   # kableExtra::add_header_above(names.head,
#   #                  bold = TRUE,
#   #                  extra_css = "padding:3px") %>%
#   kableExtra::kable_classic(full_width = TRUE,
#                 html_font = "Public Sans Light",
#                 font_size = 14) %>%
#   
#   kableExtra::column_spec(1, width = "4cm") %>%
#   kableExtra::column_spec(2:3, width = "11cm") %>%
# 
#   kableExtra::row_spec(0, bold = TRUE, extra_css = "border-bottom:1px solid black; padding:3px") %>%
#   row_spec(lines, extra_css = "bbgcolor:#d3d3d3; padding:3px") %>%
#   kableExtra::row_spec(1:base::nrow(geoda), extra_css = "padding:3px") %>%
#   kableExtra::row_spec(nrow(geoda), bold = TRUE, extra_css = "border-bottom:1px solid black; padding:3px") %>%
#   
#   # Make add_header_row cell border black
#   gsub("1px solid #ddd", "1px solid black", .)
# 
# writeLines(foo, "test.txt")

```