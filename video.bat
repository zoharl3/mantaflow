
rem r: 10, 20, 30
ffmpeg -r 30 -i frame_%%04d.png water.mp4 -y
rem pause
