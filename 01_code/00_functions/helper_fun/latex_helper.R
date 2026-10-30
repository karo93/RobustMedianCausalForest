#  add note
add_note_inside_kable <- function(kable_output, note) {
  # Pattern: matches \end{tabular} with optional }}
  pattern <- "(\\\\end\\{tabular\\}\\}*)"
  
  # Replacement: close tabular, keep any extra braces, then add vspace + minipage with note
  replacement <- paste0(
    "\\1%\n",
    "\\\\vspace{0.5em}\n",
    "\\\\begin{minipage}{0.95\\\\linewidth}\n",
    "\\\\footnotesize\n",
    "\\\\textit{Note:} ", note, "\n",
    "\\\\end{minipage}\n"
  )
  
  gsub(pattern, replacement, kable_output, perl = TRUE)
}
