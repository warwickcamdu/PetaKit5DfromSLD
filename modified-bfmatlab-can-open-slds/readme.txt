The bfmatlab here has been modified using the below procedure. The bfmatlab is version 7.3.

To open slidebook files (.sld) in MATLAB is not so easy. For some reason you can't just use the official bfmatlab (bioformats matlab toolbox). The bfmatlab here has been modified using the below procedure. 

1. Download a specific SlideBook6Reader.jar file from 3i. 3i emailed me this 28th May 2024 using the hightail file upload service. Link here, and also on CAMDU server here.

2. Install 7zip. Unzip the SlideBook6Reader.zip file. 

3. Right click the SlideBook6Reader.jar file and click 7-zip->Open archive. Move to loci\formats\in\ and drag drop the SlideBook6Reader.class onto your desktop. Move to \META-INF\lib\windows_64\ and drag drop the SBReadFile.dll onto your desktop. Go to the other directories if you have a different operating system. I have only tested Windows. 

4. Download the official bioformats matlab toolbox from here (click downloads page, then artifacts, then bfmatlab.zip). Unzip the bfmatlab.zip file. 

5. Open the new bfmatlab folder and right click on bioformats_package.jar and click 7-zip->Open archive. Move to the loci\formats\in\ and drag drop the SlideBook6Reader.class from your desktop into this location. Move to \META-INF\lib\windows_64\ and drag drop the SBReadFile.dll from your desktop into this location. 

6. Open matlab and add the bfmatlab directory to the path.

7. Open sld. files using the standard bioformats matlab method. For example: data = bfopen('C:\Users\camdu\Downloads\sld-into-matlab-working\Example-file-1T-10Z-1C.sld') 
