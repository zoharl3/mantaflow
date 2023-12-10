
@echo off
ffmpeg -r 30 -v quiet -stats -i "%~dp0\frame_%%04d.png" "%~dp0\video.mp4" -y
rem pause
