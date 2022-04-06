# NCL-script

	begin
	
	 ;============= index ===========
	  
	  var 	= z_m           		    ; replace z_m(nvar,nt) with input
									                ; nt is the temporal length of each component
	  nvar 	= dimsizes(var(:,0))
	  
    rec 	= new(nvar,"float")
	  tval 	= new(nvar,"float")
	  reg 	= new((/nvar,nyr/),"float")
	  
	  do i = 0,nvar-1
	  
	  	 rec(i)			= regline(x,var(i,:))
		 
		 b 				= rec@yintercept
		 tval(i) 		= rec@tval
		 reg(i,:)		= rec(i)*x+b
		 
	  end do 
	 

	  y 			= new((/nvar,3,nyr/),"float")
	  y(:,0,:) 		= var
	  y(:,1,:) 		= reg
	  y(:,2,:) 		= sm
  
	  ; y 			= new((/nvar,4,nyr/),"float")
	  ; y(:,0,:) = var
	  ; y(:,1,:) = reg
	  ; y(:,2,:) = sm
	  ; y(:,3,:) = clm
	  
  	  
	  tt10 		= 1.686 ; n=90
	  tt05 		= 2.024 ; n=95
	  tt01 		= 2.708 ; n=99
	  
	  tt 		= tt05
	  
	  rc 		= sqrt(tt^2/(nyr-2+tt^2))
	 
	  print("rc="+rc)
	 
	;==========plot==============
	
	
	wksname = oupath + " "
	
	wks 	= gsn_open_wks("x11",wksname)

	plot 	= new(1,graphic)


		  res                       	= True
		  res@gsnFrame             		= False
		  res@gsnDraw              		= False                 ; plot mods desired

		  res@vpXF 						= 0.15
		  res@vpYF 						= 0.6
		  res@vpWidthF  				= 0.8 
		  res@vpHeightF 				= 0.3
		  
		  res@tmXTOn                	= False
		  res@tmYROn                	= False
		  res@tmXBLabelFontHeightF 	= 0.02
		  res@tmYLLabelFontHeightF 	= 0.02
		  res@trXMinF                	= min(xtime)
		  res@trXMaxF                	= max(xtime)
		  
		  ; res@tmXBMode				= "Explicit"
		  ; res@tmXBValues				= x(1::5)
		  ; res@tmXBLabels				= xtime(1::5)
		  
		  ; res@trYMinF                	= 20
		  ; res@trYMaxF                	= 36
		  res@trYMinF                	= 24
		  res@trYMaxF                	= 33
		  
		  ; res@tiYAxisString         	= "SCSSM Intensity"
		  ; res@tiYAxisFontHeightF 		= 0.015
   
		  res@xyLineThicknesses  		= (/2.0,2.0,1.0,1.0/)
		  res@xyDashPatterns     		= (/0,1,1,1/)
		  res@xyLineColors       		= (/"black","red","blue","black"/)
 
		  res@gsnLeftStringOrthogonalPosF = 0.05
          
		  res@tiMainFontHeightF	 		= 0.025
		  res@tiMainString      	 = "SCSSM Onset Pentad"
			
          leftstring				 = (/"(a) ","(b) "/) + "RMS of 2-8-d z500 in " + (/"pre-MO","post-MO"/);"
        
		 do i = 0,nvar-1
			res@gsnLeftString      	 = leftstring(i) + " (rec="+sprintf("%5.3f",rec(i))+", tt="+sprintf("%5.3f",tval(i))+")"
			plot(i) 				 = gsn_csm_xy(wks,xtime,y(i,:,:),res)
		 end do 
			
			
			
			
			line		= True
			line@gsLineThicknessF 		= 1.5
			line@gsLineDashPattern		= 0
			line@gsLineColor			= "black"
			
			dum0 	= gsn_add_polyline(wks,plot(0),(/min(xtime),max(xtime)/),(/0,0/),line)
			dum1 	= gsn_add_polyline(wks,plot(1),(/min(xtime),max(xtime)/),(/0,0/),line)
			
	
	resp              				 	  = True
	resp@gsnPanelRowSpec 				  = True
	resp@gsnPanelCenter  				  = False
	;resp@gsnPanelTop             = 0.99
	resp@gsnPanelBottom  				  = 0.1
	resp@gsnPanelLabelBar 				= False
	resp@lbBoxLinesOn     				= False
	resp@lbLabelFontHeightF 			= 0.015
	resp@pmLabelBarHeightF 			= 0.06
	resp@gsnPanelXWhiteSpacePercent 	= 3
	resp@gsnPanelYWhiteSpacePercent 	= 3
	
	resp@txString = " "
	
	resp@txFontHeightF = 0.02	 

  gsn_panel(wks,plot,(/1,1/),resp)     

	end



