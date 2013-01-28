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

;; Transform a vector representing a row in the banmi-model for display
;;
(define (output-transform output column-config)
  (let ((is (apply append (map (lambda (type) (column-indices type column-config))
                               '(discrete ordered continuous)))))
    (map (lambda (i j) (column-inverse-transform (vector-ref output j)
                                                 (vector-ref column-config i)))
         is
         (iota (length is)))))

;; Display a header row
;;
(define (display-header column-config)
  (let* ((is (apply append (map (lambda (type) (column-indices type column-config))
                                '(discrete ordered continuous))))
         (colnames (map (lambda (i) 
                          (assq-ref (vector-ref column-config i) name:)) 
                        is)))
    (apply format #t (string-join (make-list (length is) "~a") " ") colnames)
    (newline)))

;; Returns a format string that can be used to display a row of output data
;; 
(define (output-format column-config)
  (define (n-type type)
    (length (column-indices type column-config)))
  (string-join (append (make-list (n-type 'discrete) "~d")
                       (make-list (n-type 'ordered) "~d")
                       (make-list (n-type 'continuous) "~,2f")
                       (list "~%"))
               " "))

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

;; read the configuration and allocate the model

(define argl (command-line))
(if (< (length argl) 2)
  (error "Usage: guile -s thit.scm CONFIGURATION_FILE"))

(define model-config #f)
(define data-config #f)

(do-with-file (cadr argl) read
  (lambda (form)
    (case (car form)
      ((model) (set! model-config (apply pairs (cdr form))))
      ((data) (set! data-config (primitive-eval form)))
      (else 
       (format #t "Warning: unknown section in configuration file: ~a~%" (car form))))))

(define banmi-model (new-model model-config data-config))

;; load data and perform the augmentation

(define data-file (open-input-file (assq-ref data-config file:)))

(if (assq-ref data-config header:)
  (read-line data-file))           ;; throw away the header

(do-with-file data-file read-line
  (lambda (line)
    (let* ((input (line->vector line))
           (model-input (input-transform input (assq-ref data-config columns:))))
      (apply banmi-load-row! banmi-model model-input))))

(close-input-port data-file)

(banmi-data-augmentation! banmi-model (assq-ref model-config n-iter:))

;; display the imputed data

(let* ((column-config (assq-ref data-config columns:))
       (format-string (output-format column-config))
       (imputed-data (banmi-get-imputed-data banmi-model)))
  (display-header column-config)
  (for-each (lambda (i)
              (apply format #t format-string
                     (output-transform (vector-ref imputed-data i) column-config)))
            (iota (vector-length imputed-data))))

(display "lambda ")
(for-each (lambda (x) (format #t "~,2f " x))
          (vector->list (banmi-get-lambda banmi-model)))
(newline)

(display "sigma ")
(for-each (lambda (x) (format #t "~,2f " x))
          (vector->list (banmi-get-sigma banmi-model)))
(newline)

(format #t "unique modes: ~d~%" (banmi-count-unique-modes banmi-model))
