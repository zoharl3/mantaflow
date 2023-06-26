
@echo off

if [%1]==[] goto usage

set out_path=%1

echo Copying log to %out_path%

copy _log.ans %out_path%

copy C:\prj\python\ansi2html_header.html %out_path%\_log.html
ansi2html -i < _log.ans >> %out_path%\_log.html

C:\prj\python\ansi2txt.py < _log.ans > %out_path%\_log.txt
goto :eof

:usage
@echo Usage: %0 ^<path^>
exit /B 1



