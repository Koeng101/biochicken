
(require-extension srfi-1)
(require-extension srfi-69)
(require-extension statistics)
(require-extension alist-lib)
(require-extension s)

#| 
Parameter setup
|#
(define codon_database "codon_usage.spsum")

(define codons (list "CGA" "CGC" "CGG" "CGT" "AGA" "AGG" "CTA" "CTC" "CTG" "CTT" "TTA" "TTG" "TCA" "TCC" "TCG" "TCT" "AGC" "AGT" "ACA" "ACC" "ACG" "ACT" "CCA" "CCC" "CCG" "CCT" "GCA" "GCC" "GCG" "GCT" "GGA" "GGC" "GGG" "GGT" "GTA" "GTC" "GTG" "GTT" "AAA" "AAG" "AAC" "AAT" "CAA" "CAG" "CAC" "CAT" "GAA" "GAG" "GAC" "GAT" "TAC" "TAT" "TGC" "TGT" "TTC" "TTT" "ATA" "ATC" "ATT" "ATG" "TGG" "TAA" "TAG" "TGA"))
(define amino_acids (list "R" "R" "R" "R" "R" "R" "L" "L" "L" "L" "L" "L" "S" "S" "S" "S" "S" "S" "T" "T" "T" "T" "P" "P" "P" "P" "A" "A" "A" "A" "G" "G" "G" "G" "V" "V" "V" "V" "K" "K" "N" "N" "Q" "Q" "H" "H" "E" "E" "D" "D" "Y" "Y" "C" "C" "F" "F" "I" "I" "I" "M" "W" "*" "*" "*"))

(define simple_amino_acids (delete-duplicates amino_acids))
#|
Codon table setup
|#
(define codon_list (read-lines codon_database))


(define (taxid_codons codon_list taxid)
  (if (equal? taxid (car (string-split (car codon_list) ":")))
      (cadr codon_list)
      (taxid_codons (cdr codon_list) taxid)))

(define aa_cdn (zip amino_acids codons (string-split (taxid_codons codon_list "4932"))))


(define (aa_transform aa_cdn buffer)
  (if (null? aa_cdn)
      buffer
      (let ((newlist ((lambda (data)
                  (let ((codon (cadr data))
                        (frac (cddr data)))
                    (cons codon frac))
                  )
                (car aa_cdn)))
      (aa (car (car aa_cdn)))
      (current_pair (car aa_cdn))
      )
        (if (assoc aa buffer)
            (aa_transform (cdr aa_cdn) (alist-cons aa (cons newlist (cdr (assoc aa buffer))) (alist-delete aa buffer))) ; updates buffer
            (aa_transform (cdr aa_cdn) (alist-cons aa (list newlist) buffer))) ; adds amino acid to buffer
        )))


(define codon_list (aa_transform aa_cdn '()))

(define (fctr short_list buffer)
        (if (null? short_list)
        buffer
        (fctr (cdr short_list) (+ (string->number (cadar short_list)) buffer ))))


(define (short_average short_list cdn_average buffer)
(if (null? short_list)
buffer
(short_average (cdr short_list) cdn_average (alist-update (caar short_list) (/ (string->number (cadar short_list)) cdn_average) buffer))))

(define codon_table (alist->hash-table (map (lambda (x) (cons (car x) (short_average (cdr x) (fctr (cdr x) 0) '()))) codon_list)))



(define (choose_codon table codon) (let((samples (hash-table-ref table (string codon))))
				    (random-weighted-sample 1 (alist-keys samples) (alist-values samples))))

(define (optimize_protein table amino_acids) (s-join "" (join (map (lambda (x) (choose_codon table x)) (string->list 
amino_acids)))))
