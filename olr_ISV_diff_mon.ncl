

	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
	load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

	begin
	
	inpath      = "/mnt/e/data/olr/"
	oupath      = "~/onedrive/exPrec/fig/"
  

	  latSt     = -20
	  latLt  	= 40
	  lonSt 	= 40
	  lonLt  	= 180
      yrSt  	= 1979
	  yrmd 		= 1992
	  yrLt   	= 2018     
	
	  monSt 	= 6 
	  month 	= "June"
	  
	  
	  nyr1    	= yrmd - yrSt + 1   ; 1979-1992
	  nyr2    	= yrLt - yrmd		; 1993-2018
	  
	 
	  ;============== olr ==================
	  
	  
	  ;f       	= addfile(inpath + "olr_fil.2-10d-121wgts.day.1979-2019.nc","r")
	 f       	= addfile(inpath + "olr_fil.60d-lowpass-121wgts.day.1979-2019.nc","r")
	  print(f)
	  exit
	  Time 		= f->time
	  time		= cd_calendar(Time,-5)
	  yrARR		= time(:,0)
	  monARR	= time(:,1)
	
	  lat 		= f->lat({latSt:latLt})
	  lon 		= f->lon({lonSt:lonLt})
	  nlat		= dimsizes(lat)
	  nlon		= dimsizes(lon)
	  
	  std1 		= new((/nyr1,nlat,nlon/),"float")
	  std2 		= new((/nyr2,nlat,nlon/),"float")
	  prob 		= new((/nlat,nlon/),"float")
	 
	
	  do iyr 	= yrSt,yrmd
	  
	  monId1 	:= ind(yrARR.eq.iyr .and. monARR.eq.monSt)
	  
	  Raw1 		:= f->olr_fil(monId1,{latSt:latLt},{lonSt:lonLt})
      std1(iyr-yrSt,:,:) = dim_stddev_n_Wrap(Raw1,0)
	  
	  end do 
	  
	  
	  do iyr 	= yrmd+1,yrLt
	  
	  monId2 	:= ind(yrARR.eq.iyr .and. monARR.eq.monSt)
	  
	  Raw2 		:= f->olr_fil(monId2,{latSt:latLt},{lonSt:lonLt})
      std2(iyr-yrmd-1,:,:) = dim_stddev_n_Wrap(Raw2,0)
      	  
	  end do
      
	  ave1 			= dim_avg_n_Wrap(std1,0)
	  var1 			= dim_variance_n_Wrap(std1,0)
	 
	  ave2 			= dim_avg_n_Wrap(std2,0)
	  var2 			= dim_variance_n_Wrap(std2,0)
	  
	  prob 			= ttest(ave1,var1,nyr1,ave2,var2,nyr2,False,False)
	  	  
      dt			= ave2 - ave1  
		
      copy_VarCoords(Raw1(0,:,:),dt)
	  copy_VarCoords(Raw1(0,:,:),prob)
	  
      
	 printMinMax(dt,0)
	 ;exit
	;==========plot==============
	
	
    wksname 	= oupath+"olr-60d-low_diff_"+month+"_95"
	
	wks 		= gsn_open_wks("x11",wksname)
	
	gsn_define_colormap(wks,"precip_diff_12lev")
	
	plot 		= new(1,graphic)


	res          				= True
	res@gsnDraw 				= False
	res@gsnFrame 				= False
	res@gsnAddCyclic 			= False

	res@cnInfoLabelOn         	= False
	res@cnFillOn              	= True               ; turn on color for contours
	res@cnLinesOn            	= False             ; turn off contour lines
	res@cnLineLabelsOn       	= False              ; turn off contour line labels
    
	; lonSt 		= -100
	; lonLt 		= 180
	res@mpCenterLonF 			= (lonSt+lonLt)/2.
	res@mpMaxLatF 				= latLt
	res@mpMinLatF		 		= latSt
	res@mpMinLonF 				= lonSt
	res@mpMaxLonF 				= lonLt

	res@mpGeophysicalLineColor   	   = "black"
	res@mpGeophysicalLineThicknessF = 1.5

	res@tmXTOn    				= False
	res@tmYROn    				= False
	res@tmXBLabelFontHeightF 	= 0.02
	res@tmYLLabelFontHeightF 	= 0.02
	res@tmXBTickSpacingF  		= 30

	res@gsnSpreadColors 		= True

	res@pmLabelBarOrthogonalPosF  = 0.15
	res@lbLabelFontHeightF       	= 0.02
	res@lbLabelBarOn             	= True
	res@lbBoxLinesOn             	= False
	;res@lbOrientation				= "vertical"

	res@cnLevelSelectionMode 	= "ExplicitLevels"
	res@cnFillColors      		= (/2,3,4,5,6,7,9,10,11,12,13,14/)
	
	res@gsnLeftStringFontHeightF = 0.022
	res@gsnRightString 			= " "
    
    nn 	= (/"(a)","(b)","(c)","(d)"/)    ;"
 	
	res@cnLevels  			= (/-5,-4,-3,-2,-1,0,1,2,3,4,5/)*1.0

	  
	   res@gsnLeftString 	= "(a) Diff. olr_60d-low in "+ month+" (1993-2018 minus 1979-1992)"
	   plot         	= gsn_csm_contour_map_ce(wks,dt,res)
	   
	 

	delete(res@mpCenterLonF)
	delete(res@mpMaxLatF)
	delete(res@mpMinLatF)
	delete(res@mpMaxLonF)
	delete(res@mpMinLonF)
	delete(res@mpGeophysicalLineColor)
	delete(res@mpGeophysicalLineThicknessF)
	delete(res@cnFillColors )
	
	
	delete(res@cnLevels)
	res@cnFillOn 						= False
	res@cnLinesOn           			= False              ; turn off contour lines
	res@cnLineLabelsOn       			= False
	res@cnLevels 				       	= (/0.05,0.1/)
	res@cnFillDotSizeF 				= 0.003
	res@cnFillColor 					= "gray40"
	res@gsnLeftString					= " "
	
	shade                  				= True
	shade@gsnShadeFillType 				= "pattern"
	shade@gsnShadeLow 	   				= 17
	;shade@gsnShadeHigh 					= 17
	
	
	
	   plot0         = gsn_csm_contour(wks,prob,res)
	   plot0         = gsn_contour_shade(plot0,0.05,0.1,shade)
	   overlay(plot,plot0)
     
	   
	 
  
	   
	line 					= True
	line@gsLineColor 		= "purple"
	line@gsLineThicknessF	= 2.0
	
	latS 		= 5 
	latN		= 22
	lonW		= 105
	lonE		= 130
	xpts 		= (/lonW,lonE,lonE,lonW,lonW/)
	ypts		= (/latS,latS,latN,latN,latS/)
	
	; dum 		= gsn_add_polyline(wks,plot(0),xpts,ypts,line)
	; dum0 		= gsn_add_polyline(wks,plot(1),xpts,ypts,line)
	
	draw(plot)
	frame(wks)
	
	end


