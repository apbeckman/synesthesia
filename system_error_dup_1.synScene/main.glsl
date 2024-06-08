vec4 renderMain(void) {
	vec2 position = _uv;
	position += _muv;
	vec3 col = vec3(0.);
	vec4 newcol = vec4(1.);
	
	
    	vec2 pxl = 1./RENDERSIZE;
    	vec2 screenwarp = sin(_uv.x*10. + TIME*1.)*pxl*2.*syn_Presence;
    	screenwarp += _noise(_uvc*600. + 4.*sin(TIME))*pxl*10.;
    	screenwarp.x += sin(screenwarp.x*20.)/RENDERSIZE.x*1.;
    	screenwarp.y += 10.*pxl.y*_noise(_pixelate(sin(_uv*.5 + TIME*1.3),10.) * 40. + TIME*0.2)*syn_Level;
    	
    	
    	float vert_glitch = sin(_uv.y*20.+syn_BassTime) * sin(_uv.y*30.+syn_MidTime) * sin(_uv.y*35.+syn_Time) * sin(_uv.x+_uv.y+TIME) * sin(_uv.x*2.+_uv.y*10.+TIME);
        vert_glitch = smoothstep(0.,0.7,vert_glitch);
        
        vec2 sqr_uv = _rotate(_uvc,0.5 + sin(TIME*0.1));
        float sqr_glitch = step(0.5,sin(sqr_uv.x * 8.));
        sqr_glitch += step(0.5,sin(sqr_uv.y * 8.+syn_BPMTwitcher));
        sqr_glitch = mod(sqr_glitch,1.99);
        sqr_uv = _rotate(sqr_uv,TIME*0.1);
        sqr_glitch += step(0.5,sin(sqr_uv.x * 3.+syn_Time));
        sqr_glitch = mod(sqr_glitch,1.99);
        sqr_uv = _rotate(sqr_uv,TIME*0.1);
        sqr_glitch += step(0.3,sin(sqr_uv.y * 6. + TIME));
        sqr_glitch = mod(sqr_glitch,1.99);
        
        screenwarp += 5. * vert_glitch / RENDERSIZE;
        screenwarp += 5. * sqr_glitch / RENDERSIZE;
        
        
    	screenwarp *= twitchy*3.;
    	if (warp_beat > 0.5 && syn_BPMConfidence > 0.35) screenwarp *= syn_BPMTri2;
    	
    	screenwarp.x = _pixelate(screenwarp.x,pxl.x);
    	screenwarp.y = _pixelate(screenwarp.y,pxl.y);
    	
    	screenwarp.x += _uvc.y*(syn_BPMSin2*0.5 - syn_BPMConfidence/2.)*0.25 * pow(window_boogie, 2.);
    	screenwarp.y += _uvc.x*(syn_BPMSin4*0.5 - syn_BPMConfidence/2.)*0.25 * pow(window_boogie, 2.);
    	
    	screenwarp += _uv;
    
    if (PASSINDEX == 0) {
        
	    float mouse_dist = 0.;
    	vec2 lastuv = _uv;
    	vec2 pxl = 1./RENDERSIZE;
	    lastuv = screenwarp;
    	if (warp > 0.2 || (warp_beat > 0.5 && syn_BPMConfidence > 0.5)) {
            lastuv.x += (_noise(_uvc * 40. * warp_noise_scale + TIME) * 0.01 - 0.005) * warp_noise;
            lastuv.y += (_noise(_uvc * 40. * warp_noise_scale + TIME - 1000.) * 0.01 - 0.005) * warp_noise;
            vec2 hitswarp = vec2(
                syn_BassHits * 0.01 * syn_BassLevel * syn_Level - 0.001 * smoothstep(0.2,0.24,sin(_uv.x * 100. + syn_Time * 200.)) - 0.001,
                syn_HighHits * 0.01 * syn_HighLevel * syn_Level - 0.001 * smoothstep(0.2,0.24,sin(_uv.x * 100. + syn_Time * 200.)) - 0.001
            );
            lastuv += hitswarp * warp_hits;
        	lastuv.x = floor(lastuv.x*RENDERSIZE.x)/RENDERSIZE.x;
        	lastuv.y = floor(lastuv.y*RENDERSIZE.y)/RENDERSIZE.y;
    	}
	    
	    vec2 glitch_amt = vec2(0.);
	    for (float i=0.; i<3.; i++) {
	        glitch_amt.x += _sqPulse(_uv.y, sin(syn_Time*(1.+i)/20.)*0.5 + 0.5, fract(PI*TIME*0.05)*0.1);
	        glitch_amt.y += _sqPulse(_uv.x, sin(syn_HighTime*(1.+i)/80.)*0.5 + 0.5, fract(PI*TIME*0.2)*0.1);
	    }
        if (syn_ToggleOnBeat > 0.5) glitch_amt.x *= -1.;
        if (syn_RandomOnBeat > 0.5) glitch_amt.y *= -1.;
	    lastuv.x += _pixelate(glitch_amt.x*glitch_amt.x*constant_warp*0.1, pxl.x);
	    lastuv.y += _pixelate(glitch_amt.x*glitch_amt.x*constant_warp*0.1, pxl.y);
    	vec3 lastcol = texture(passA,lastuv).rgb;
        lastcol = texelFetch(passA,ivec2(lastuv*RENDERSIZE),0).rgb;
	
    	if (warp > 0.2) {
	        lastcol = _rgb2hsv(lastcol);
	        lastcol.r += sin(TIME);
    	    lastcol = _hsv2rgb(lastcol);
	
	        lastuv = _rotate(lastuv, TIME*0.1);
	
	        lastuv.x += (lastcol.r - lastcol.g) * 0.01 * syn_BassHits * warp_melt;
	        lastuv.xy += (lastcol.r - lastcol.g) * 0.015 * syn_MidHighHits * warp_melt;
	
	        lastuv = _rotate(lastuv, TIME*-0.1);
	        
	        lastuv = fract(lastuv);
	        ivec2 lpuv = ivec2(lastuv*RENDERSIZE);
	        lastcol = texelFetch(passA,lpuv,0).rgb;
	    }
	    col = lastcol;
        if (warp > 0.5) {
            col = _rgb2hsv(col);
            // col.g += 0.002;
            if (platform > 0.5 && platform < 1.5)
                col.g = 0.;
            // col.b -= 0.003 * lastcol.g;
            col = _hsv2rgb(col);
        }
	
    	if (crash > 0.5) {
            if (platform < 0.5) { // win xp   
    	        col = texture(winxp_bsod, screenwarp).rgb;
            } else if (platform < 1.5) { // mac os7
    	        col = texture(os7_crash, screenwarp).rgb;
            }
    	}
	
    	if (welcome > 0.5 || FRAMECOUNT < 3) {
            if (platform < 0.5) { // win xp   
    	        col = texture(winxp_welcome, screenwarp).rgb;
            } else if (platform < 1.5) { // mac os7
    	        col = texture(os7_happy, screenwarp).rgb;
            }
    	}
	
    	// draw error windows
    	if (_click.x > 0.5) {
            vec2 drawuv = (_uvc - _muvc) - window_size;
    	    drawuv *= 1. / window_size;
    	    drawuv -= 0.5;
    	    if (drawuv.x < -1. && drawuv.x > -2. && drawuv.y < -1. && drawuv.y > -2.) {
                if (platform < 0.5) { // win xp   
    	            newcol = texture(winxp_error, drawuv);
                    col = mix(col,newcol.rgb,newcol.a);
                } else if (platform < 1.5) { // mac os7
                    drawuv.y *= 3.;
                    if (drawuv.y < -4. && drawuv.y > -5.) {
    	                newcol = texture(os7_error, drawuv);
                        col = mix(col,newcol.rgb,newcol.a);
                    }
                }
    	    }
    	}
	
	    // random error windows
    	if (error_blast > 0.1 || error_once > 0.5) {
            vec2 drawuv = (_uvc - vec2(_rand(TIME)*2.-1.,_rand(100.+TIME)*1.5-.75)) - window_size;
    	    drawuv *= 1. / window_size;
    	    drawuv -= 0.5;
    	    if (drawuv.x < -1. && drawuv.x > -2. && drawuv.y < -1. && drawuv.y > -2.) {
                if (platform < 0.5) { // win xp   
    	            newcol = texture(winxp_error, drawuv);
                    col = mix(col,newcol.rgb,newcol.a);
                } else if (platform < 1.5) { // mac os7
                    drawuv.y *= 3.;
                    if (drawuv.y < -4. && drawuv.y > -5.) {
    	                newcol = texture(os7_error, drawuv);
                        col = mix(col,newcol.rgb,newcol.a);
                    }
                }
    	    }
    	}
	
	    // stack of error windows
	    if (error_stack > 0.2 && mod(FRAMECOUNT,6)<1) {
	        vec2 offset = vec2(TIME)*0.3;
	        vec2 ar = vec2(RENDERSIZE.x/RENDERSIZE.y, RENDERSIZE.y/RENDERSIZE.x);
	        offset = mod(offset, 2.*vec2(1.,ar.y))-vec2(1.,ar.y);
            vec2 drawuv = _uvc + offset +- window_size;
	        drawuv *= 1. / window_size;
	        drawuv -= 0.5;
	        if (drawuv.x < -1. && drawuv.x > -2. && drawuv.y < -1. && drawuv.y > -2.) {
                if (platform < 0.5) { // win xp   
	                newcol = texture(winxp_error, drawuv);
                    col = mix(col,newcol.rgb,newcol.a);
                } else if (platform < 1.5) { // mac os7
                    drawuv.y *= 3.;
                    if (drawuv.y < -4. && drawuv.y > -5.) {
                        newcol = texture(os7_error, drawuv);
                        col = mix(col,newcol.rgb,newcol.a);
                    }
                }
            }
    	}
    }
    
    if (PASSINDEX==1) {
        vec3 a = texture(passA,screenwarp).rgb;
        // vec3 lp = texture(syn_FinalPass,_uv).rgb;
        // col = mix(a,lp,0.9);
        
        col = a;
        
        if (media_player > 0.5) {
            vec2 wuv  = screenwarp; // popup window uv
            wuv.y += 0.03;
            // wuv.y += sin(TIME*6.)*0.03*(0.5+syn_Level);
            wuv.y += mix(sin(TIME*6.),syn_BPMSin2,syn_BPMConfidence)*0.03*(0.5+syn_Level);
            wuv -= 0.5;
            wuv *= 2.2;
            wuv += 0.5;
            wuv = mix(wuv, _uv*0.8 + 0.1, fullscreen_media);
            if (platform < 0.5) { // win xp
                if (wuv.x<0.903 && wuv.x>0.097 && wuv.y < 0.955 && wuv.y > 0.093) {
                    col = vec3(0.1373, 0.3373, 0.8745)*0.8;
                }
                if (wuv.x<0.9 && wuv.x>0.1 && wuv.y < 0.95 && wuv.y > 0.1) {
                    col = vec3(0.1373, 0.3373, 0.8745);
                    // top right button outline
                    if(wuv.y<0.946 && wuv.y > 0.904 && wuv.x > 0.867 && wuv.x < 0.896) col = vec3(1.);
                    // top right button
                    if(wuv.y<0.945 && wuv.y > 0.905 && wuv.x > 0.868 && wuv.x < 0.895) col = vec3(0.902, 0.4235, 0.3137);
                    // todo: derive sdf of X for close button
                    vec2 close_btn_uv = wuv;
                    _uv2uvc(close_btn_uv);
                    close_btn_uv = abs(close_btn_uv);
                    float close_btn_x = length(close_btn_uv-min(close_btn_uv.x+close_btn_uv.y,0.015)*.5) + 0.01;
                        
                    col = mix(col,vec3(1.),smoothstep(0.013,0.012,close_btn_x));
                }
                if (wuv.x<0.9 && wuv.x>0.1 && wuv.y < 0.9 && wuv.y > 0.1) {
                    if (_exists(syn_UserImage)) {
                        col = _textureUserImage((wuv - 0.5)*1.25 + vec2(0.5,0.5)).rgb; 
                    } else {
                        col = texture(syn_FinalPass,(wuv - 0.5)*1.1 + 0.5).rgb;
                    }
                }
            } else if (platform < 1.5) { // os7
                if (wuv.x<0.903 && wuv.x>0.097 && wuv.y < 0.955 && wuv.y > 0.093) {
                    col = vec3(.3);
                }
                if (wuv.x<0.9 && wuv.x>0.1 && wuv.y < 0.95 && wuv.y > 0.1) {
                    col = vec3(.85);
                    // top stripes
                    if(wuv.y<0.945 && wuv.y > 0.905) col = vec3(smoothstep(0.,0.05,sin(wuv.y*920.))) * 0.6 + 0.4;
                    // top left button
                    if(wuv.y<0.945 && wuv.y > 0.905 && wuv.x > 0.13 && wuv.x < 0.154) col = vec3(.9);
                    // top right button
                    if(wuv.y<0.945 && wuv.y > 0.905 && wuv.x > 0.846 && wuv.x < 0.87) col = vec3(.9);
                    // top left button inset
                    if(wuv.y<0.943 && wuv.y > 0.904 && wuv.x > 0.132 && wuv.x < 0.153) col = vec3(0.6784, 0.6706, 0.7843) * 0.7;
                    if(wuv.y<0.942 && wuv.y > 0.908 && wuv.x > 0.133 && wuv.x < 0.15) col = vec3(0.6784, 0.6706, 0.7843);
                    // top right button inet
                    if(wuv.y<0.943 && wuv.y > 0.904 && wuv.x > 0.848 && wuv.x < 0.868) col = vec3(0.6784, 0.6706, 0.7843) * 0.7;
                    if(wuv.y<0.942 && wuv.y > 0.908 && wuv.x > 0.849 && wuv.x < 0.863) col = vec3(0.6784, 0.6706, 0.7843);
                    
                }
                if (wuv.x<0.9 && wuv.x>0.1 && wuv.y < 0.9 && wuv.y > 0.1) {
                    if (_exists(syn_UserImage)) {
                        col = texture(syn_UserImage,(wuv - 0.5)*1.25 + vec2(0.5,0.5)).rgb; 
                    } else {
                        col = texture(syn_FinalPass,(_uv - 0.5)*1.1 + 0.5).rgb;
                    }
                }
            }
        }
        
        if (platform > 0.5 && platform < 1.5) { // os7 grayscale
            col = col.ggg;
        }
    }
	return vec4(col, 1.);
}