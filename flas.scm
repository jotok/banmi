(load-extension "./libthit" "banmi_thit")

(define *max-rows* 300)
(define *bds-discrete* '(2 2 2 2 3))
(define *n-continuous* 4)
(define *dp-weight* 50)
(define *lambda-a* 2)
(define *lambda-b* 8)

(define model (new-model *max-rows* *bds-discrete* *n-continuous*
                         *dp-weight* *lambda-a* *lambda-b*))

(display (get-lambda model))
(newline)
