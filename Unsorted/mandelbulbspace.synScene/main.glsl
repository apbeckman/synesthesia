bool flop = (flop_bool > 0.5); 
bool flip = (flip_bool > 0.5); 


////////////////////////////////////////////////////////////
// Equirec_z33d+  by mojovideotech
//
// spherical mod of :
// interactiveshaderformat.com/\2129
// based on :
// shadertoy.com/\XlXcWj
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////

#define 	twpi  	6.283185307179586  	// two pi, 2*pi
#define 	pi   	3.141592653589793 	// pi

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);
  
	float k = 0.0, T = smoothTime*0.1;
	vec2 R = RENDERSIZE.xy;
	float TW = (tan(smoothTime/(3.5+warpfastorslow))/2.5);
	float IF = (((smoothTimeC)/10.)/5.);
  	vec2 M = vec2((mXY.x+(infinity2*IF)+(timewarp*TW)),mXY.y+infinity*IF)*R.xy;
  	vec2 sph = (_xy.xy / RENDERSIZE.xy) * vec2(twpi, pi);
    vec3 rd = vec3(sin(sph.y) * sin(sph.x), cos(sph.y), sin(sph.y) * cos(sph.x)); 
    vec3 uv;
    if (flip) {
  			uv = rd.zyx;
  			if (flop) uv = rd.zxy;
  		} 
  		else if (flop) uv = rd.yzx;
  		else uv = rd.xyz;
  	for (float i = 0.0; i < 24.0; i++) {
    	if (i > depth) { break; }
    	vec3 p = vec3(uv.xy + center.xy, k-1.0);
    	float a = T;
    	p.zy *= mat2(cos(rXY.y+a), -sin(rXY.y+a), sin(rXY.y+a), cos(rXY.y+a));
	 	a /= 2.0;
    	p.yx *= mat2(cos(rXY.x+a), -sin(rXY.x+a), sin(rXY.x+a), cos(rXY.x+a));
    	a /= 2.0;
    	p.zx *= mat2(cos(rZ+a), -sin(rZ+a), sin(rZ+a), cos(rZ+a));
    	vec3 z = p;
    	float c = 2.0;
    	for (float j = 0.0; j < 9.0; j++) {
      		float r = length(z);
        	if (r > 6.0) { 
        		k += log(r) * r / c / pi;
				break;
        	}
      		float d = acos(z.z / r) * (zoom + (zoom+6.0) * M.x / R.x);
      		float b = atan(z.y, z.x) * (zoom + (zoom+6.0)  * M.y / R.y);
      		c = pow(r, 7.0) * 5.0 * c / r + 1.0;
      		z = pow(r, 7.0) * vec3(sin(d) * cos(b), -sin(d) * sin(b), -cos(d)) + p;
    	}
    	out_FragColor = vec4(1.0 - i / (3.0*glow) - k + p / (12.0-saturation), 1.0);
      	if (log(length(z)) * length(z) / c < e) {
       		break;
    	}
  	}

return out_FragColor; 
 } 
