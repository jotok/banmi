(use-modules (ice-9 rdelim)
             (ice-9 format)
             (srfi srfi-11))

(load-extension "./libthit" "banmi_thit")

(define *max-rows* 300)
(define *n-iter* 50)

(define *bds-discrete* '(2 2 2 2 3))
(define *n-continuous* 4)

(define *dp-weight* 50)
(define *lambda-a* 2)
(define *lambda-b* 8)

;; utility functions

(define (line->numbers line)
  (map string->number
       (delete "" (string-split line #\space))))

(define (do-lines fn)
  (do ((line (read-line) (read-line)))
      ((eof-object? line))
      (fn line)))

(define (continuous->ordered x max-value)
  (inexact->exact (floor (* (1+ max-value) x))))

(define (ordered->continuous o max-value)
  (/ (+ o 0.5) (1+ max-value)))

;; run the multiple imputation

(define model (new-model *max-rows* *bds-discrete* *n-continuous*
                         *dp-weight* *lambda-a* *lambda-b*))

;; This function will be called to load rows from flas1.txt, a space-delimited
;; file with 14 columns of data. See flas_description.txt for more information.
(define (load-data-from-file)
  (read-line) ;; throw away the header
  (do-lines 
    (lambda (line)
      (let-values (((lan2 lan3 lan4 age pri sex mlat flas satv satm eng hgpa cgpa grd)
                    (apply values (line->numbers line))))
        (load-row! model 
                   lan2 lan3 lan4 
                   (if (>= age 0) (1- age) -1) 
                   (if (>= pri 0) (1- pri) -1) 
                   (if (>= satv 0) (ordered->continuous satv 800) -1)
                   (if (>= satm 0) (ordered->continuous satm 800) -1)
                   (if (>= hgpa 0) (/ hgpa 4.0) -1) 
                   (if (>= cgpa 0) (/ cgpa 4.0) -1))))))

(with-input-from-file "data/flas1.txt" load-data-from-file)
(data-augmentation! model *n-iter*)

;; show the result
(for-each (lambda (row)
            (format #t "~2d ~2d ~2d ~2d ~2d ~4d ~4d ~4,2f ~4,2f~%"
                    (vector-ref row 0)
                    (vector-ref row 1)
                    (vector-ref row 2)
                    (1+ (vector-ref row 3))
                    (1+ (vector-ref row 4))
                    (continuous->ordered (vector-ref row 5) 800)
                    (continuous->ordered (vector-ref row 6) 800)
                    (* (vector-ref row 7) 4)
                    (* (vector-ref row 8) 4)))
          (vector->list (get-imputed-data model)))

(display "lambda ")
(for-each (lambda (x) (format #t "~,2f " x))
          (vector->list (get-lambda model)))
(newline)

(display "sigma ")
(for-each (lambda (x) (format #t "~,2f " x))
          (vector->list (get-sigma model)))
(newline)

(format #t "unique modes: ~d~%" (count-unique-modes model))
