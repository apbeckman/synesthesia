float line_segment(in vec2 p, in vec2 a, in vec2 b) {

	vec2 ba = b - a;

	vec2 pa = p - a;

	float h = clamp(dot(pa, ba) / dot(ba, ba), 0., 1.);

	return length(pa - h * ba);

}



vec4 renderMain(void) {

	vec2 position = _uv;

	position += _mouse.xy / RENDERSIZE;

	vec3 col = vec3(0.);

	float mouse_dist = 0.;

	vec3 newcol = vec3(1.);

	float line = 1.;

	if (_click.x > 0.5) {
        float brush_sizer = brush_size*brush_size*syn_Level;
    	line = smoothstep(
    	    brush_sizer-brush_smooth*brush_size,
    	    brush_sizer+brush_smooth*brush_size,
    	    line_segment(_uvc, _muvc, lastmouse)
    	  );
    	newcol = vec3(line);

	}

	vec2 lastuv = _uv;

    lastuv += _noise(_uvc * 20. + TIME) * 0.0006 - 0.0003;

	vec3 lastcol = texture(syn_FinalPass,lastuv).rgb;


	// col = mix(col * 5., lastcol, 0.9997 * line);
	col = lastcol;
    col = _rgb2hsv(col);
    col.r += 0.005;
    col.g += 0.001;
    col.b -= 0.001;
    col = _hsv2rgb(col);
	
	if (_click.x > 0.5) {
	    
        if (draw_mode < 0.5) {
            
            col = mix(col, vec3(0.3,0.5,0.8), 1. - 0.9995 * line);
            
        } else if (draw_mode < 1.5) {
            
            vec2 drawuv = (_uv - _muv) * vec2(1.,-1.) - brush_size;

	        drawuv *= 1. / brush_size;
	        
	        drawuv -= 0.5;
	        
	        vec3 newcol = texture(syn_UserImage, drawuv).rgb;

	        vec3 newcolhsv = _rgb2hsv(newcol);
	        vec3 colhsv = _rgb2hsv(col);
	        if (newcolhsv.b > 0.1)
    	        if (colhsv.b < newcolhsv.b)
            	   col = mix(col,newcol,(1.-line)*(newcol.r+newcol.g+newcol.b)/3.);
        }
	}

	return vec4(col, 1.);

}

