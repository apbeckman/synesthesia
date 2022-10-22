

			//******** Common Code Begins ********



#define pal(a,b,c,d,e) ((a) +(b)*sin((c)*(d) + (e)))

#define rot(a) mat2(cos(a), -sin(a), sin(a), cos(a))


// TRIANGLE MODE - substitutes sine for tri 
//#define sin(a) (asin(sin(a))*1.)


mat3 getOrthogonalBasis(vec3 direction){
    direction = normalize(direction);
    vec3 right = normalize(cross(vec3(0,1,0),direction));
    vec3 up = normalize(cross(direction, right));
    return mat3(right,up,direction);
}
float cyclicNoise(vec3 p, bool turbulent, float time){
    float noise = 0.;
    
    float amp = 1.;
    float gain = 0.8 + sin(p.z*0.2)*0.2;
    const float lacunarity = 1.4;
    const int octaves = 5;
    
    const float warp = 1.2;    
    float warpTrk = 1.5 ;
    const float warpTrkGain = .2;
    
    vec3 seed = vec3(-4,-2.,0.5);
    mat3 rotMatrix = getOrthogonalBasis(seed);
    
    for(int i = 0; i < octaves; i++){
        
        p += sin(p.zxy*warpTrk + vec3(0,-time*2.,0) - 2.*warpTrk)*warp; 
        noise += sin(dot(cos(p), sin(p.zxy + vec3(0,time*0.3,0))))*amp;
    
        p *= rotMatrix;
        p *= lacunarity;
        
        warpTrk *= warpTrkGain;
        amp *= gain;
    }
    
    if(turbulent){
        return 1. - abs(noise)*0.5;
    
    }{
        return (noise*0.25 + 0.5);

    }
}


float plaIntersect( in vec3 ro, in vec3 rd, in vec4 p )
{
    return -(dot(ro,p.xyz)+p.w)/dot(rd,p.xyz);
}



// Cyclic Noise from nimitz. 
// fwidth analytic SDF AA suggested by Fabrice (another solution suggested before that by mla, who afaik found it from iq)

// check out triangle mode in common


const float slices = 124.;

const float sliceDepth = 2.;

float yOffs = 1.5;

const float width = 0.;

const float aaSteps = 1.; // aa unused


float fun(vec3 p){
    float f = 0.;
    
    p.z += TIME;
    f += cyclicNoise(p*2. + vec3(0,TIME*0.2,0), true, TIME*0.1)
        - mix(-0.5,2.,smoothstep(0.,1.,-p.y))
        + mix(2.,0.,smoothstep(0.,1.,-abs(p.y )+ yOffs ));
        
    f = max(f,-dot(p.xy,p.xy) + 0.);


    return f;
}




vec3 get(in vec2 fragCoord){
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    vec3 col = vec3(0);
    
    vec3 ro = vec3(0,yOffs,0);
    
    vec3 rd = normalize(vec3(uv,3));
    rd.yz *= rot(0.4);
    
    for(float i = 0.; i < slices; i++ ){
        
        float t = plaIntersect( ro - vec3(0,i/slices*sliceDepth,0.), rd, vec4(0,1,0,1.) );
        t = abs(t);
        vec3 p = ro + rd*t;
        
        float fn = fun(p);
        
        
        float fw = fwidth(fn) ;
        
        
        vec3 c = pal(0.5,vec3(0.6,0.05,.15),vec3(4.9,1.,9),1.,i/slices*22. + smoothTimeC*0.5);
        
        float d = (fn)/fw;
        
        float g = smoothstep(1.,0.,abs(d)- width);
        float gb = smoothstep(1.,0.,d - width);
        
        
        col = mix(col, c*c*0.7*1. + 0.*smoothstep(1.51,0.,abs(fn)*20.0005), gb);
        col = mix(col, c*2., g);
        col = max(col,0.);


         
        
        /*
            col = mix(col, c*c*0.7, gb);
            col = mix(col, c*2., g);
        
        }*/
        
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
    
    
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    //col = mix(col,col*col*0.1,exp(-dot(uv,uv)*5422.));
    
    col = max(col, 0.);
	
    
    col = pow(col,vec3(0.4545));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}