
@echo off

set out_path=c:\prj-external-libs\mantaflow\out\

copy _log.ans %out_path%

copy C:\prj\python\ansi2html_header.html %out_path%\_log.html
ansi2html -i < _log.ans >> %out_path%\_log.html

C:\prj\python\ansi2txt.py < _log.ans > %out_path%\_log.txt

