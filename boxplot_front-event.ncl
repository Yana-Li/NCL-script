
  
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

	begin


	inpath      	= "~/onedrive/exPrec/data/"
	oupath      	= "~/onedrive/exPrec/fig/"
      
	
      yrSt  		= 1979
	  yrmd 			= 1992
	  yrLt   		= 2018 
	  
      nyr 			= yrLt - yrSt + 1	  
	  nyr1    		= yrmd - yrSt + 1
	  nyr2    		= yrLt - yrmd
	  nd 			= 30
      
	  nt 			= 5
	  Raw 			= new((/nyr*nd,nt/),"float")
	  ;Raw  			= asciiread(inpath + "frontIndex/front-event_daily_Jun_1979-2018_20-35N_105-125E_ERA5.txt",(/nyr*nd,nt/),"float")
	  Raw  			= asciiread(inpath + "frontIndex/front-event-weight_daily_Jun_1979-2018_20-35N_105-125E_ERA5.txt",(/nyr*nd,nt/),"float")
	  Raw@_FillValue 	= -999
	  
	  nevn 			= 10
	  evnNum 		= new((/nyr/),"float")
	  evnDay1		= new((/nyr1,nevn/),"float",-999)
	  evnDay2		= new((/nyr2,nevn/),"float",-999)
	  speed1		= new((/nyr1,nevn/),"float",-999)
	  speed2		= new((/nyr2,nevn/),"float",-999)
	  
	  do i 			= 0,nyr-1	  
	     evnNum(i)  	= max(Raw(i*nd:(i+1)*nd-1,4)) 			 
	  end do 
	  
	  do i 			= 0,nyr1-1
	     do j 		= 1,evnNum(i)
		    ind0 				:= ind(Raw(i*nd:(i+1)*nd-1,4).eq.j)
			day0 				:= Raw(i*nd:(i+1)*nd-1,1)
			speed0 				:= Raw(i*nd:(i+1)*nd-1,2)
			
			evnDay1(i,j-1) 		= day0(ind0(0))
	        speed1(i,j-1) 		= speed0(ind0(0))
			
	     end do 
	  end do 
	  
	  do i 			= nyr1,nyr-1
	     do j 		= 1,evnNum(i)
		 
		    ind0 				:= ind(Raw(i*nd:(i+1)*nd-1,4).eq.j)
			day0 				:= Raw(i*nd:(i+1)*nd-1,1)
			speed0 				:= Raw(i*nd:(i+1)*nd-1,2)
			
			evnDay2(i-nyr1,j-1) 			= day0(ind0(0))
	        speed2(i-nyr1,j-1) 			= speed0(ind0(0))
			
	     end do 
	  end do 
	  
	  printMinMax(evnDay1,-1)
	  printMinMax(speed1,-1)
	  
	  
      evnDay1_1D 			= ndtooned(evnDay1)
	  evnDay2_1D 			= ndtooned(evnDay2)
	  speed1_1D 			= ndtooned(speed1)
	  speed2_1D 			= ndtooned(speed2)
	  
	  evnDay1_1D@_FillValue = -999
	  evnDay2_1D@_FillValue = -999
	  speed1_1D@_FillValue 	= -999
	  speed2_1D@_FillValue 	= -999
	  
	  var1 					= evnDay1_1D
	  var2 					= evnDay2_1D
	  
	  indx     				= ind(.not.ismissing(var1))
	  varx     				= var1(indx)
	  indy     				= ind(.not.ismissing(var2))
	  vary     				= var2(indy)
	  
	  ;---------- box values -----------
	  
	  y 					= new((/2,5/),"float")
	  mean  				= new((/2/),"float")
	  mean(0)  				= dim_avg(varx)
	  mean(1)  				= dim_avg(vary)
	  
	  i = 0
       x1 := varx
       nx  = dimsizes(x1)
	   print(nx)
       ix := dim_pqsort(x1,2)
	   
       y(i,0)	=  x1(toint(floor(nx*0.1))-1) 
       y(i,1)	=  x1(toint(floor(nx*0.25))-1)
       y(i,2)	=  dim_median(x1)
       y(i,3)	=  x1(toint(floor(nx*0.75))-1) 
       y(i,4)	=  x1(toint(floor(nx*0.9))-1) 
     	 
       i = 1
       x1 := vary
       nx  = dimsizes(x1)
	   print(nx)
       ix := dim_pqsort(x1,2)
	   
       y(i,0)	=  x1(toint(floor(nx*0.1))-1) 
       y(i,1)	=  x1(toint(floor(nx*0.25))-1)
       y(i,2)	=  dim_median(x1)
       y(i,3)	=  x1(toint(floor(nx*0.75))-1) 
       y(i,4)	=  x1(toint(floor(nx*0.9))-1) 
     	 
    	   
	   
	   
	
	wksName 		= oupath+"boxplot_front-event_days"
	wks 			= gsn_open_wks("x11",wksName)       ; open file to plot 
	plot 			= new(1,graphic)
   
 
	res                      	= True
	res@gsnFrame             	= False
	res@gsnDraw             	= False                 ; plot mods desired
	res@vpWidthF 				= 0.5
	res@vpHeightF 				= 0.6
 
	res@tmXTOn             		= False
	res@tmYROn             		= False
	res@tmXBLabelFontHeightF 	=	0.02
	res@tmYLLabelFontHeightF 	=	0.02
	res@tiYAxisFontHeightF    =   0.02
	res@tiXAxisFontHeightF    =   0.02
	res@tiYAxisString 			=   "Days"
    res@tmXBLabels 				= (/"1979-1992","1993-2018"/)
    res@gsnLeftStringFontHeightF = 0.04
	
	; res@trYMaxF 				= 2.0
	; res@trYMinF 				= 0.0
	
	llres                   	= True			
    llres@gsLineThicknessF  	= 3.0
	
	opti          				= True			
    opti@boxWidth 				= 0.3		          ; width in unit of xaxis		
    opti@boxColors 				= (/"blue","red"/)  	
	
	x  = (/0,1/)
	res@tiMainString   		= "Days of Front Events in June"
    plot = boxplot(wks,x,y,opti,res,llres)
	
	mres               = True                     ; marker mods desired
    mres@gsMarkerIndex = 1                        ; polymarker style
    mres@gsMarkerSizeF = 0.02                      ; polymarker size
	mres@gsMarkerThicknessF = 2.0
    mres@gsMarkerColor = "black"                  
  

   dum1 = gsn_add_polymarker(wks,plot,x,mean,mres)
   
   ; mres@gsMarkerIndex = 3
   ; mres@gsMarkerColor = "black"   
   ; dum2 = gsn_add_polymarker(wks,plot,x,pre_2020,mres)
	
	draw(plot)
	frame(wks)

end
