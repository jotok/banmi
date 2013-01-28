(model
  max-rows: 300
  n-iter: 50
  
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
           (sex ignore)
           (mlat ignore)
           (flas ignore)
           (satv ordered max: 800)
           (satm ordered max: 800)
           (eng ignore)
           (hgpa continuous max: 4)
           (cgpa continuous max: 4)
           (grd ignore))
