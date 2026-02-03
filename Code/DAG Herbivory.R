library(dagitty)
# library(ggdag)

Caterpillar = dagitty("dag{
  temperature -> caterpillar;
  precipitation -> caterpillar;
  Plant_L -> caterpillar;
  temperature -> Plant_L;
  precipitation -> Plant_L;
  
    Plant_L [unobserved]
}"
)

plot(Caterpillar)

# (minimum) adjustment sets.
adjustmentSets(Caterpillar, exposure = "temperature", outcome = "caterpillar")
 