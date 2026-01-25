
# Plant genera selection ----
# Use "All" if you want to return all values in the set

get_genera = function(plantGenus, fullDataset) {
  if ("All" %in% plantGenus) {
    unique(fullDataset$plantGenus)
  } else {
    plantGenus
  }
}


# Site selection function----
get_site = function(site, fullDataset) {
  if ("All" %in% site) {
    unique(fullDataset$Name)
  } else {
    site
  }
}

