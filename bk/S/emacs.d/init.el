;; .emacs

;;; uncomment this line to disable loading of "default.el" at startup
;; (setq inhibit-default-init t)

;; enable visual feedback on selections
;(setq transient-mark-mode t)

;; default to better frame titles
(setq frame-title-format
      (concat  "%b - emacs@" (system-name)))

;; default to unified diffs
(setq diff-switches "-u")

;; always end a file with a newline
;(setq require-final-newline 'query)

;;; uncomment for CJK utf-8 support for non-Asian users
;; (require 'un-define)

(set-background-color "#333366")

(set-foreground-color "#ffffff")

(set-cursor-color "#ffffff")

(global-font-lock-mode t)

(setq font-lock-maximum-decoration t)
(setq fast-lock nil)
(setq lazy-lock nil)
(setq jit-lock t)

(tool-bar-mode 0)

(set-scroll-bar-mode 'right)

(global-set-key "\C-h" (quote delete-backward-char))

(if (eq window-system '(x))
    (progn
      (global-set-key [end] 'end-of-buffer )
      (global-set-key [home] 'beginning-of-buffer )
    )
)

(setq display-time-day-and-date t)
(display-time)

(setq fill-column 80)
(setq text-mode-hook 'turn-on-auto-fill)
(setq default-major-mode 'text-mode)

(setq inhibit-startup-message t)

(setq make-backup-files nil)

(put 'upcase-region 'disabled nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; KUMAC mode -- get it from http://uther1.phy.ornl.gov/~morrison/
;;
(autoload 'kumac-mode "/home/exp/Analys/users/ss/emacs.d/kumac-mode.el" "Mode for KUMAC files." t)
(setq auto-mode-alist (cons '("\\.kumac$" . kumac-mode) auto-mode-alist))
(add-hook 'kumac-mode-hook (lambda () (abbrev-mode 1)))
(add-hook 'kumac-mode-hook (lambda () (font-lock-mode 1)))
(setq font-lock-maximum-decoration t)
(global-font-lock-mode t)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ANALYS mode
;;
(autoload 'ana-mode "/home/exp/Analys/users/ss/emacs.d/ana-mode.el" "Mode for editing ANA files." t)
(setq auto-mode-alist (cons '("\\.ana$" . ana-mode) auto-mode-alist))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; ANALYS PRM mode
;;
(autoload 'prm-mode "/home/exp/Analys/users/ss/emacs.d/prm-mode.el" "Mode for editing ANA files." t)
(setq auto-mode-alist (cons '("\\.prm$" . prm-mode) auto-mode-alist))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; INC files --> FORTRAN Mode
;;
(setq auto-mode-alist (cons '("\\.inc$" . fortran-mode) auto-mode-alist))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; fh files --> FORTRAN Mode
;;
(setq auto-mode-alist (cons '("\\.fh$" . fortran-mode) auto-mode-alist))