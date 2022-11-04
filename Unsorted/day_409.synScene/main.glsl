vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// Fork of "day idk wip" by jeyko. https://shadertoy.com/view/ttKcWc
// 2021-01-29 16:17:36

// Fork of "Day 407" by jeyko. https://shadertoy.com/view/3tycWV
// 2021-01-29 08:24:51

// Fork of "Day 405" by jeyko. https://shadertoy.com/view/ttVcRt
// 2021-01-27 11:17:18
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.25;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.25;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.25;


const float slices = 100.;

const float aaSteps = 1.; // not really steps, it's the exponentially ^3 area growing area around the fragCoord 

const float disp = .4;

const float width = 0.002;

#define pal(a,b,c,d,e) ((a) +(b)*sin((c)*(d) + (e)))

#define sin(x) mix( sin(x), abs(fract((x + 3.14/2.*sin(smoothTimeC*0.5))/3.14*1.) - 0.5)*2., (sin(smoothTimeC)*0.5 + 1.)*0.52)

float fun(float p, float py){
    py *= 90. + sin(smoothTime)*15.;
    
    py += smoothTime*8.;
    p *= 1.;
    float f = abs(sin(p*9. + sin(py*(.15)*1.5  )*9.));
    
    f = pow(max(f,0.001),3.);
    
    f += (sin(py*0.1 + smoothTimeC + sin(p*8. + smoothTimeC*.5 + sin(py*2.)*0.1)))*.2;
    
    return f*disp;
}


float graph(float y, float fn0){
  return smoothstep(1. ,0., 
                    abs(fn0-y)/fwidth(fn0-y)- width);
}
float graphNoAbs(float y, float fn0){
  return smoothstep(1.,0., 
                    -(fn0-y)/fwidth(fn0-y) - width);
}


vec3 get(in vec2 fragCoord){
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    vec3 col = vec3(0);
    
    vec2 ouv = uv;
    uv = vec2(atan(uv.y,uv.x),log(length(uv) ));
    
    for(float i = 0.; i < slices; i++ ){
        vec2 p = uv + vec2(0.,i/slices*5. );
        
        float fn = fun(p.x,i/slices);
        
        float bg = graphNoAbs( p.y, fn);
        
        col *= 1.- bg;
        
        vec3 c = pal(0.5,vec3(0.6,0.05,.15),vec3(4.9,1.,9),1.,i + smoothTime*8.);
        col += 2.*bg*max(c,0.);
        
        c -= 1.;
        
        #define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
        
        c.xz *= rot(-1.5);
        
        c += 1.;
        col = mix(col, c*0.2*c*c, graph( p.y, fn ));
        
    }
    
    
    return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 col = vec3(0);
    
    
    for(float i =0.; i < aaSteps*aaSteps + min(float(FRAMECOUNT),0.)   ; i++){
    	col += get(fragCoord + vec2(mod(i,aaSteps),floor(i/aaSteps))/aaSteps);
    }
    col /= aaSteps*aaSteps;
    
    
    col = max(col, 0.);
	
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    col = mix(col,col*col*0.1,exp(-dot(uv,uv)*5422.));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}