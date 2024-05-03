PRO movie_example

outputdir = '/scratch/sepid/DATA/2017.08.12/postprocess/aligned/'
videodir = '/home/seki2695/OUTPUT/video/'

  Compile_Opt idl2

     ; Open the video recorder.
     video_file = videodir+'sample2.avi'
     video = IDLffVideoWrite(video_file, Format='avi')
     
     ; Configure the video output for the PNG files you
     ; plan to add.
     framerate = 5
     
     xsize = 272*5
     ysize = 183*5
     stream = video.AddVideoStream(xsize, ysize, framerate)
     
     ; Get the files.
     files = File_Search(outputdir+'crispex.stokes.6302*.png', COUNT=fileCnt)
     stop
     FOR j=0,fileCnt-1 DO BEGIN
        
       ; Read the PNG file.
       image = Read_PNG(files[j])
;stop
       ; Save the image in the video stream
       void = video.Put(stream, image)
  ;stop   
     ENDFOR
     
     ; Clean things up.
     video.Cleanup
stop
END
