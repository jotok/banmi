(use-modules (ice-9 format)
             (srfi srfi-11))
(use-syntax (ice-9 syncase))

(read-set! keywords 'postfix)

;; Make an alist from a vararg.
(define (pairs . args)
  (let loop ((acc '()) (args args))
    (if (< (length args) 2)
      (begin (if (= (length args) 1)
               (format #t "Warning: unpaired value: ~a~%" (car args)))
             (reverse acc))
      (let-values (((k v) (apply values (list-head args 2))))
        (loop (cons (cons k v) acc)
              (list-tail args 2))))))

(define (vector-multi-ref vec indices)
  (let loop ((acc '()) (is indices))
    (if (null? is)
      (reverse acc)
      (loop (cons (vector-ref vec (car is)) acc) (cdr is)))))

;; Transform the data section of the configuration file to an alist.
(define-syntax data
  (syntax-rules ()
    ((_ columns: (name type . options) ...)
     (list (cons columns: (vector (pairs name: 'name type: 'type . options) ...))))
    ((_ k v . rest)
     (cons (cons k v) (data . rest)))))

(define (column-indices type column-config)
  (filter (lambda (i) (eq? (assq-ref (vector-ref column-config i) type:) 
                           type))
          (iota (vector-length column-config))))

;; Open the file and apply fn to each value returned by read-fn.
(define (do-with-file file read-fn fn)
  (with-input-from-file file
    (lambda ()
      (do ((datum (read-fn) (read-fn)))
          ((eof-object? datum))
          (fn datum)))))

;; Load the confiuration and determine the banmi configuration

(define model-config #f)
(define data-config #f)

(do-with-file "config.scm" read
  (lambda (form)
    (case (car form)
      ((model) (set! model-config (apply pairs (cdr form))))
      ((data) (set! data-config (primitive-eval form)))
      (else 
       (format #t "Warning: unknown section in configuration file: ~a~%" (car form))))))
