(use-modules (ice-9 format)
             (ice-9 rdelim)
             (srfi srfi-11))
(use-syntax (ice-9 syncase))

(load-extension "./libthit" "banmi_thit")

(read-set! keywords 'postfix)

;; Make an alist from a vararg.
;;
(define (pairs . args)
  (let loop ((acc '()) (args args))
    (if (< (length args) 2)
      (begin (if (= (length args) 1)
               (format #t "Warning: unpaired value: ~a~%" (car args)))
             (reverse acc))
      (let-values (((k v) (apply values (list-head args 2))))
        (loop (cons (cons k v) acc)
              (list-tail args 2))))))

;; Takes a list of vector indices and returns a list of the corresponding
;; values at those indices.
;;
(define (vector-multi-ref vec indices)
  (let loop ((acc '()) (is indices))
    (if (null? is)
      (reverse acc)
      (loop (cons (vector-ref vec (car is)) acc) (cdr is)))))

;; Merge default options into an alist
;;
(define (merge-defaults alist defaults)
  (if (null? defaults)
    alist
    (merge-defaults (if (assq-ref alist (caar defaults))
                      alist
                      (cons (car defaults) alist))
                    (cdr defaults))))

;; Transform a column description to an alist
;;
(define-syntax column-description-transform
  (syntax-rules()
    ((_ (name type . options))
     (merge-defaults (pairs name: 'name type: 'type . options)
                     '((min: . 0))))))

;; Transform the data section of the configuration file to an alist.
;;
(define-syntax data
  (syntax-rules ()
    ((_ columns: desc ...)
     (list (cons columns: (vector (column-description-transform desc) ...))))
    ((_ k v . rest)
     (cons (cons k v) (data . rest)))))

;; Return a list of column indices corresponding to columns of the given type.
;;
(define (column-indices type column-config)
  (filter (lambda (i) (eq? (assq-ref (vector-ref column-config i) type:) 
                           type))
          (iota (vector-length column-config))))

;; Given a column description, return a function that transforms an input datum
;; to the form expected by the model.
;;
(define (column-transform x column)
  (let* ((min (assq-ref column min:))
         (size (- (assq-ref column max:) min)))
    (case (assq-ref column type:)
      ((discrete)
       (- x min))
      ((ordered)
       (/ (+ (- x min) 0.5) (1+ size)))
      ((continuous)
       (/ (- x min) size)))))

;; Same as column-transform, but negative values pass through as -1.
;;
(define (column-transform-1 x column)
  (if (< x 0)
    -1
    (column-transform x column)))

;; Given a column description, return a function that transforms a data
;; column back to its original value range
;;
(define (column-inverse-transform x column)
  (let* ((min (assq-ref column min:))
         (size (- (assq-ref column max:) min)))
    (case (assq-ref column type:)
      ((discrete)
       (+ x min))
      ((ordered)
       (+ (inexact->exact (floor (* x (1+ size)))) min))
      ((continuous)
       (+ (* x size) min)))))

;; Allocate a new banmi model with the given configuration
;;
(define (new-model model-config data-config)
  (let* ((column-config (assq-ref data-config columns:))
         (discrete-columns (column-indices 'discrete column-config))
         (ordered-columns (column-indices 'ordered column-config))
         (continuous-columns (column-indices 'continuous column-config))
         (n-continuous (+ (length ordered-columns) (length continuous-columns)))
         (bds-discrete (map (lambda (col) (1+ (- (assq-ref col max:)
                                                 (assq-ref col min:))))
                            (vector-multi-ref column-config discrete-columns))))
    (new-banmi-model (assq-ref model-config max-rows:)
                     bds-discrete
                     n-continuous
                     (assq-ref model-config dp-weight:)
                     (assq-ref model-config lambda-a:)
                     (assq-ref model-config lambda-b:))))

;; Transform a vector of input values to a list of values that can be loaded
;; into the model
;;
(define (input-transform input column-config)
  (define (input-transform-type type)
    (map (lambda (i) (column-transform-1 (vector-ref input i) 
                                         (vector-ref column-config i)))
         (column-indices type column-config)))
  (apply append 
         (map input-transform-type '(discrete ordered continuous))))

;; Apply fn to each line read from file.
;;
(define (do-with-file file read-fn fn)
  (define (proc port)
    (do ((datum (read-fn port) (read-fn port)))
        ((eof-object? datum))
        (fn datum)))
  (if (input-port? file)
    (proc file)
    (call-with-input-file file proc)))

;; Convert a space-delimited string to a list of numbers
;;
(define (line->vector line)
  (list->vector 
   (map string->number (delete "" (string-split line #\space)))))

;;
;; Main program
;;

(define model-config #f)
(define data-config #f)

(do-with-file "config.scm" read
  (lambda (form)
    (case (car form)
      ((model) (set! model-config (apply pairs (cdr form))))
      ((data) (set! data-config (primitive-eval form)))
      (else 
       (format #t "Warning: unknown section in configuration file: ~a~%" (car form))))))

(define banmi-model (new-model model-config data-config))
