
//#define mx (10.*iMouse.x/iResolution.y)
//#define iTime (iTime + 2.4 + mx)
vec3 getRd(vec3 ro,vec3 lookAt,vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0), dir));
	vec3 up = normalize(cross(dir, right));
	return dir + right*uv.x + up*uv.y;
}
#define dmin(a,b) a.x < b.x ? a : b
#define pmod(p,x) mod(p, x) - x*0.5

vec4 r14c(float i){return texture(image30, vec2(i));}

float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;    
    return max(p.x, max(p.y, p.z));
}
#define pi acos(-1.)
#define tau (2.*pi)
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define tunnRotRate

vec2 id;
vec2 map(vec3 p){
	vec2 d = (vec2(10e7));

    p.xy *= rot(0. + p.z*0.1 + 0.1*TIME);
    
    
    for (float i = 0.; i < 4.; i++){
    	p = abs(p);
    	p.xy *= rot((inOrOut*(0.4+Angle))*pi);
        p.x -= 0.2;
        p.x *= 1. + 0.4*atan(p.x, p.y)/pi;
        //p.y += 0.1;
    }    

    p.xy -= 2.0;
    
    
    p.y = abs(p.y);
    
    
    p.y -= 1. + sin(TIME*0.1)*0.2;
    
    #define modSz 0.75
    id = floor(p.xz/modSz);
    //vec2 
    
    p.xy -= 0.8;
    p.xz = pmod(p.xz, modSz);
    
    for (float i = 0.; i < (5); i++){ //separation, rangge 1-6
    	p = abs(p);
        //  p.y -= 0.28 - sin(TIME*0.2)*0.08- 0.1;
        p.x += 0.04;
    	p.xy *= rot(0.6*pi + id.y*6.  + 0.9);
        if (i == 3.5){
        	p.xz *= rot((script_time+(script_high_time*0.35)) + id.y); //rotate the thingies
        }
    }     

    d = dmin(d, vec2(sdBox(p, vec3(modSz*0.3 + sin(TIME*0.26)*0.1)), 1.)); 
    
    d.x *= 0.5;
    return d;
} 
/*
vec3 getNormal(vec3 p){
	vec2 t = vec2(0.001,0);
    return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}*/
    
vec3 glow = vec3(0);

//#define pal(q,w,e,r,t) (q + w*sin(e*r + t))
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 q = fragCoord.xy / RENDERSIZE.xy;
	vec2 uv = (-RENDERSIZE.xy + 2.0*fragCoord.xy)/RENDERSIZE.y;
	
/* void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0);
*/
    vec3 col = vec3(0);
    vec3 ro = vec3(0.,0,0);
    ro.z += (smooth_basstime*0.25+script_time)*1.25;
    
    float rate = ro.z*0.1+0.1*TIME;
    
    ro.xy += vec2(sin(rate), cos(rate))*2.;
    
    vec3 lookAt = ro + vec3(0,0,4);
    float rotRate = script_time*0.05 + sin(TIME*0.13);
    lookAt.xz += vec2(
    	sin(rotRate),
    	cos(rotRate)
    );
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    vec3 p = ro; float t = 0.;
    for (int i = 0; i < 165; i++){
    	vec2 d = map(p);
								#define pal(q,w,e,r,t) (q + w*cos( tau*(e*r + t))
        //glow += exp(-d.x*70.)* pal(vec3(0.5,0.6,0.7)*1., 0.35, id.y*0.2 + iTime*0.4 + 1.*p.z*(sin(iTime)*0.001), vec3(0.4, 0.9,0.2), 0. + p.z*0.02));
        //glow += exp(-d.x*20.)* pal(vec3(0.5,0.6,0.7)*1., 0.45, id.y*0.2 + iTime*0.4 + 0.*p.z*(sin(iTime)*0.001), vec3(0.4, 0.9,0.2), 0. + p.z*0.02));
        //zglow += exp(-d.x*20.)* pal(vec3(0.5,0.6,0.7)*0.2, 0.95, id.y*0.05 + iTime*1. + 0.*p.z*(sin(iTime)*0.001), vec3(0.1, 0.9,0.2), 0.5 + p.z*0.02));
        //glow += exp(-d.x*60.)* pal(0.5, 0.45, id.y*0.2 + iTime*2., vec3(0.1, 0.4,0.8), 0.5)) ;
        glow += exp(-d.x*40.)* pal(1., 0.35, id.y*0.01 + (script_high_time*0.125+(TIME*0.175))*.1, vec3(0.4, 0.5,0.9), 0.9 + p.z*0.012)) ;
        if(d.x < 0.00025){
            /*
            vec3 n = getNormal(p);
            vec3 l = normalize(vec3(1));
            vec3 h = normalize(l - rd);
            float diff = max(dot(n,l),0.);
            float spec = max(dot(n,h),0.);
            float fres = pow(1. - max(dot(n,-rd), 0.),5.);
            */
            //col += fres*diff*3.;
            
        	break;
        }
        if (t > 85.){
        	break;
        }
        t += d.x;
        p = ro + rd*t;
    }	
    
    //float bass = pow(texture(iChannel1, vec2(0.,0.14)).x, 4.);
    
    col += glow*(0.01);
    
    col = mix(col, vec3(0), pow(clamp(t*.002 - 0.1, 0., 1.), 2.));
    col = smoothstep(0.,1., col);
    //col = smoothstep(0.,1., col);
    
    //col = pow(col , vec3(1.8,1.0,1.));
    
    
    //col.g = pow(col.g, 2. - 0.5*( col.r + col.b*0.1));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}