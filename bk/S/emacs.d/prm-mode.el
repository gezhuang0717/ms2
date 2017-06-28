;;; prm-mode-el -- Major mode for Prmlys code

;; Author: Gibelin Julien <gibelin@rarfaxp.riken.go.jp>
;; Created: May 2003
;; Keywords: PRM major-mode

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
(defvar prm-mode-hook nil)
(defvar prm-mode-map nil
  "Keymap for PRM major mode.")

(if prm-mode-map nil
  (setq prm-mode-map (make-keymap)))


(defconst prm-font-lock-keywords-1
  (list
   ; These define the beginning and end of each PRM entity definition
   '("^c.*" 0 font-lock-string-face)
   '("^0\n" 0 font-lock-warning-face)
   '("^1\n" 0 font-lock-warning-face)
   '("\\," 0 font-lock-warning-face))
  "Minimal highlighting expressions for PRM mode.")

(defvar prm-font-lock-keywords prm-font-lock-keywords-1
  "Default highlighting expressions for PRM mode.")


(defvar prm-mode-syntax-table nil
  "Syntax table for prm-mode.")

(defun prm-create-syntax-table ()
  (if prm-mode-syntax-table
	  ()
	(setq prm-mode-syntax-table (make-syntax-table))
	
    ; This is added so entity names with underscores can be more easily parsed
	(modify-syntax-entry ?_ "w" prm-mode-syntax-table)
  
	; Comment styles are same as C++
	(modify-syntax-entry ?/ ". 124b" prm-mode-syntax-table)
	(modify-syntax-entry ?* ". 23" prm-mode-syntax-table)
	(modify-syntax-entry ?\n "> b" prm-mode-syntax-table))
  (set-syntax-table prm-mode-syntax-table))

(defun prm-mode ()
  "Major PRM files."
  (interactive)
  (kill-all-local-variables)
  (prm-create-syntax-table)
  
  ;; Set up font-lock
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults
		'(prm-font-lock-keywords))
  
  ;; Register our indentation function
;;  (make-local-variable 'indent-line-function)
;;  (setq indent-line-function 'prm-indent-line)
  
  (setq major-mode 'prm-mode)
  (setq mode-name "PRM")
  (run-hooks 'prm-mode-hook))

(provide 'prm-mode)

;;; prm-mode.el ends here



