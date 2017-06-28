;;; ana-mode-el -- Major mode for Analys code

;; Author: Gibelin Julien <gibelin@rarfaxp.riken.go.jp>
;; Created: May 2003
;; Keywords: ANA major-mode

;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License as
;; published by the Free Software Foundation; either version 2 of
;; the License, or (at your option) any later version.

;; This program is distributed in the hope that it will be
;; useful, but WITHOUT ANY WARRANTY; without even the implied
;; warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
;; PURPOSE.  See the GNU General Public License for more details.

;; You should have received a copy of the GNU General Public
;; License along with this program; if not, write to the Free
;; Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
;; MA 02111-1307 USA

;;; Code:
(defvar ana-mode-hook nil)
(defvar ana-mode-map nil
  "Keymap for ANA major mode.")

(if ana-mode-map nil
  (setq ana-mode-map (make-keymap)))


(defconst ana-font-lock-keywords-1
  (list
   ; These define the beginning and end of each ANA entity definition
   ; "analys" "gate" "or" "and" "xygate" "hst1" "hst2" "c" "profile"
   '("^analys\\>" 0 font-lock-keyword-face )
   '("^\\(and\\|gate\\|or\\|xygate\\|leff\\|stop\\|hst.?\\|e.*\\)" 1 font-lock-keyword-face ) 
   '("^\\(prof\\|profile\\|profi\\|profil\\)\\b" 1 font-lock-keyword-face ) 
   '("^c.*" 0 font-lock-comment-face)
   '("\\," 0 font-lock-warning-face)
   '("'.*'" 0 font-lock-string-face))
  "Minimal highlighting expressions for ANA mode.")

(defvar ana-font-lock-keywords ana-font-lock-keywords-1
  "Default highlighting expressions for ANA mode.")


(defvar ana-mode-syntax-table nil
  "Syntax table for ana-mode.")

(defun ana-create-syntax-table ()
  (if ana-mode-syntax-table
	  ()
	(setq ana-mode-syntax-table (make-syntax-table))
	
    ; This is added so entity names with underscores can be more easily parsed
	(modify-syntax-entry ?_ "w" ana-mode-syntax-table)
  
	; Comment styles are same as C++
	(modify-syntax-entry ?/ ". 124b" ana-mode-syntax-table)
	(modify-syntax-entry ?* ". 23" ana-mode-syntax-table)
	(modify-syntax-entry ?\n "> b" ana-mode-syntax-table))
  (set-syntax-table ana-mode-syntax-table))

(defun ana-mode ()
  "Major Analys files."
  (interactive)
  (kill-all-local-variables)
  (ana-create-syntax-table)
  
  ;; Set up font-lock
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults
		'(ana-font-lock-keywords))
  
  ;; Register our indentation function
;;  (make-local-variable 'indent-line-function)
;;  (setq indent-line-function 'ana-indent-line)
  
  (setq major-mode 'ana-mode)
  (setq mode-name "ANA")
  (run-hooks 'ana-mode-hook))

(provide 'ana-mode)

;;; ana-mode.el ends here



