(model
  max-rows: 300
  n-iter: 250
  n-imputations: 10
  
  dp-weight: 50
  lambda-a: 2
  lambda-b: 8)

(data 
  file: "data/flas1.txt"
  header: #t
  columns: (lan2 discrete max: 1)
           (lan3 discrete max: 1)
           (lan4 discrete max: 1)
           (age discrete min: 1 max: 2)
           (pri discrete min: 1 max: 3)
           (sex discrete max: 1)
           (mlat ordered max: 110)
           (flas ordered max: 40)
           (satv ordered max: 800)
           (satm ordered max: 800)
           (eng ordered max: 120)
           (hgpa continuous max: 4)
           (cgpa continuous max: 4)
           (grd discrete min: 1 max: 2))
