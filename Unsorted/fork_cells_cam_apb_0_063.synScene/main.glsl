

			//******** BuffA Code Begins ********

#define hue(h) clamp( abs( fract(h + vec4(2,1,4,0)/1.) * 6. - 3.) -1. , 0., 1.)

vec2 rand( vec2 p ) {
    return fract(sin(vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3))))*43758.5453);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
  
    //stuff to expose
    float size = 1.5*(1.0+Size);
    const float uvFac = 10.;
    const float colFac = .5;
    
    
    vec2 ouv = fragCoord/RENDERSIZE.xy;        
    vec2 uv = (fragCoord - RENDERSIZE.xy*.5) / (RENDERSIZE.y*size);    
    vec2 luv = uv;
    
    vec4 texIn = texture(syn_UserImage, ouv);
    vec2 mp = texIn.rb;
    
    uv *= 100. + sin(smoothTimeC*.5+mp.x*uvFac);
   
    vec2 iuv = floor(uv);
    vec2 guv = fract(uv);      

    float mDist = 10.;
   
    vec3 col = vec3(.1);
       
    for (float y= -1.; y <= 1.; y++) {
        for (float x= -1.; x <= 1.; x++) {            
            vec2 neighbor = vec2(x, y);            
            vec2 point = rand(iuv + neighbor);
            point = .5 + .5*sin(smoothTimeC*2. + 6.2831*point);
            vec2 diff = neighbor + point - guv;            
            float dist = length(diff);                      
           
            mDist = min(mDist, dist);                        
        }
    } 
       
    float l = length(luv);    
    col = hue(fract(mDist*.95 + smoothTimeB*.1 + l + mp.x*colFac)).rgb;
    fragColor = vec4(col,1.0)*.05 + texture(BuffA, ouv) *.95;    
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = fragCoord/RENDERSIZE.xy;

    fragColor = texture(BuffA, uv);
        
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}