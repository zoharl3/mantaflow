
@echo off
rem r: 10, 20, 30, 50
ffmpeg -v quiet -stats -r 50 -i "%~dp0\frame_%%04d.png" "%~dp0\water.mp4" -y
rem pause
