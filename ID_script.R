# Characters allowed in IDs
chars <- c(
  LETTERS[!LETTERS %in% c("O", "I")],
  as.character(2:9)
)

# Generate one random 5-character ID
make_id <- function() {
  paste0(sample(chars, 5, replace = TRUE), collapse = "")
}

# Generate 5000 unique IDs
set.seed(123)  # optional: makes the result reproducible

ids <- character(0)

while (length(ids) < 5000) {
  new_ids <- replicate(5000 - length(ids), make_id())
  ids <- unique(c(ids, new_ids))
}

# Put them in a data frame
id_list <- data.frame(
  individual_id = ids
)

# Check
nrow(id_list)                    # should be 5000
length(unique(id_list$individual_id))  # should also be 5000

# Save to CSV
write.csv(
  id_list,
  "individual_IDs.csv",
  row.names = FALSE
)