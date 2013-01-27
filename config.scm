(model
  max-rows: 300
  n-iter: 50
  
  dp-weight: 50
  lambda-a: 2
  lambda-b: 8)

(data 
  file: "data/flas1.txt"
  header: #t
  columns: (lan2 discrete)
           (lan3 discrete)
           (lan4 discrete)
           (age discrete min: 1)
           (pri discrete min: 1)
           (sex ignore)
           (mlat ignore)
           (flas ignore)
           (satv ordered max: 800)
           (satm ordered max: 800)
           (eng ignore)
           (hgpa continuous max: 4)
           (cgpa continuous max: 4)
           (grd ignore))
