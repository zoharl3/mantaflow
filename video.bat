
@echo off
rem r: 10, 20, 30
ffmpeg -r 30 -i %~dp0\frame_%%04d.png %~dp0\water.mp4 -y
rem pause
