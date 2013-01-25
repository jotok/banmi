(load-extension "./libthit" "banmi_thit")

(define *max-rows* 300)
(define *bds-discrete* '(2 2 2 2 3))
(define *n-continuous* 4)
(define *dp-weight* 50)
(define *lambda-a* 2)
(define *lambda-b* 8)

(define model (new-model *max-rows* *bds-discrete* *n-continuous*
                         *dp-weight* *lambda-a* *lambda-b*))

(load-row! model 1 2 3 4 5 1.2 3.4 5.6 7.8)

(display (get-imputed-data model))
(newline)
