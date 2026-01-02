generaGoodSpatial = c("Acer", "Aesculus", "Amelanchier", "Betula", "Carpinus", "Carya", 
            "Celtis", "Cercis", "Cornus", "Crataegus", "Diospyros", "Elaeagnus", 
            "Fagus", "Fraxinus", "Hamamelis", "Ilex", "Juglans", "Kalmia", 
            "Lindera", "Liquidambar", "Liriodendron", "Lonicera", "Magnolia", 
            "Malus", "Morella", "Morus", "Nyssa", "Platanus", "Prunus", 
            "Quercus", "Rhododendron", "Sassafras", "Tilia", "Ulmus", 
            "Vaccinium", "Viburnum")



min.Genus.nSurv = 5
julianWindow = 120:230
min.nJulianWeekYearSite = 7
min.nSurvWeekYearSite = 20
min.nYear = 3
TempDayWindow = 120: 230
minLatRange = 5

goodGenusSurv = fullDataset %>% 
  filter( julianweek %in% julianWindow,
    !HerbivoryScore %in% c(-128, -1),
    plantGenus %in% generaGoodSpatial
  ) %>% 
  group_by(Name, Latitude, Longitude, Year, ObservationMethod, julianweek, plantGenus) %>% 
  summarise(nSurv.genus = n_distinct(ID)) %>% 
  filter(nSurv.genus >= min.Genus.nSurv) %>% 
  as.data.frame()


goodGenusYear = goodGenusSurv %>% 
  group_by(Name, Latitude, Longitude, ObservationMethod, Year, plantGenus) %>% 
  summarise(nJulianWeek = n_distinct(julianweek)) %>% 
  filter(nJulianWeek >= min.nJulianWeekYearSite) %>% 
  as.data.frame()


goodGenusYear %>% 
  group_by(plantGenus) %>% 
  summarise(latMin = min(Latitude),
            latMax = max(Latitude),
            latRange = round((latMax - latMin), 2)) %>% 
  filter(latRange >= minLatRange) # A dosen good plant!
