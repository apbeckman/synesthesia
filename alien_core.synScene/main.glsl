//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// GLKITTY 2016.

vec3 rotateY(vec3 v, float t){
    float cost = cos(t); float sint = sin(t);
    return vec3(v.x * cost + v.z * sint, v.y, -v.x * sint + v.z * cost);
}

float smin( float a, float b, float k )
{
    float h = clamp(( 0.5+test)+(0.5+test)*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float noise(vec3 p){
    
    float t = smoothTime;
    vec3 np = normalize(p);
    
    // kind of bi-planar mapping
    float a = texture(image48,(t/20.+TIME*0.00125)+np.xy).x;      
    float b = texture(image48,(t/20.+TIME*0.00125)+np.yz).x;
    
    a = mix(a,.5,abs(np.x));
    b = mix(b,.5,abs(np.z));
    
    float noise = a+b-.4;    
    noise = mix(noise,.5,abs(np.y)/2.);
        
    return noise;
}

float map(vec3 p){
    
    // spheres
    float d = (-1.*length(p)+3.)+1.5*noise(p);    
    d = min(d, (length(p)-1.5)+1.5*noise(p) );  
    
    // links
    float m = 1.5; float s = .03;    
    d = smin(d, max( abs(p.x)-s, abs(p.y+p.z*.2)-.07 ) , m);          
    d = smin(d, max( abs(p.z)-s, abs(p.x+p.y/2.)-.07 ), m );    
    d = smin(d, max( abs(p.z-p.y*.4)-s, abs(p.x-p.y*.2)-.07 ), m );    
    d = smin(d, max( abs(p.z*.2-p.y)-s, abs(p.x+p.z)-.07 ), m );    
    d = smin(d, max( abs(p.z*-.2+p.y)-s, abs(-p.x+p.z)-.07 ), m );
    
    return d;
}

float color( vec3 p){
   return 0.; 
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    
    // Ray from UV
	vec2 uv = fragCoord.xy * 2.0 / RENDERSIZE.xy - 1.0;
    uv.x *= RENDERSIZE.x / RENDERSIZE.y;
    vec3 ray = normalize(vec3(1.*uv.x,1.*uv.y,1.));
    
    // Color    
    vec3 color = vec3(0);    
    const int rayCount = 2000;
    
    // Raymarching
    float t = 0.;
    for (int r = 1; r <= rayCount; r++)
    {
        // Ray Position
        vec3 p = vec3(0,0,-3.) + ray * t;        
        
        // Rotation 
       	p = rotateY(p, rotatexy.x/RENDERSIZE.x * 2.* 3.14);  
        p = rotateY(p,((smoothTime*0.125))/3.);
        
        // Deformation 
    	float mask = max(0.,(1.-length(p/3.)));
    	p = rotateY(p,mask*sin((smoothTimeC)/10.)*1.2);        
        p.y += sin(smoothTime*0.1+p.x)*mask*.015;
        p *= 1.1+(sin((smoothTimeC*0.0125))*mask*.63);

        // distance
        float d =  map(p);   
        
        //color
        if(d<0.001 || r == rayCount )
        {                 
            
            float iter = float(r) / float(rayCount);
            float ao = (1.-iter);
            ao*=ao;
            ao = 1.-ao;
                        
            float mask = max(0.,(1.-length(p/2.)));            
            mask *= abs(sin(smoothTime*-1.5+length(p)+p.x)-.2);            
            color += 2.*vec3(.1,1.,.8) * max(0.,(noise(p)*4.-2.6)) * mask;            
            color += vec3(.1,.5,.6) * ao * 6.;            
            color += vec3(.27,.2,.4)*(t/8.);
                       
            color *= 2.;
            color -= .15;
                        
            break;          
        }
        
        // march along ray
        t +=  d *.5;        
    }
    
    // vignetting effect by Ippokratis
    // https://www.shadertoy.com/view/lsKSWR
    uv = fragCoord.xy / RENDERSIZE.xy;
    uv *=  1.0 - uv.yx; 
    float vig = uv.x*uv.y * 20.0;    
    vig = pow(vig, 0.25);        
    color *= vig;
    
    //color adjustement
    color.y *= .8;
    color.x *= 1.5;
    
	fragColor = vec4(color, 1);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}