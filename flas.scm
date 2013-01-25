(use-modules (ice-9 rdelim)
             (srfi srfi-11))

(load-extension "./libthit" "banmi_thit")

(define *max-rows* 300)
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
        (load-row! model lan2 lan3 lan4 age pri satv satm hgpa cgpa)))))

(with-input-from-file "data/flas1.txt" load-data-from-file)

(display (get-data model))
