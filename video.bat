
@echo off
rem r: 10, 20, 30, 50, 60
ffmpeg -r 50 -v quiet -stats -i "%~dp0\frame_%%04d.png" "%~dp0\water.mp4" -y
rem pause
