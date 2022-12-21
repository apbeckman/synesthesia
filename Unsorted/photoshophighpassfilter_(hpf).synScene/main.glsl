

// Fork of "Photoshop Highpass Filter (HPF)" by mbeytekin. https://shadertoy.com/view/tldcRn
// 2022-12-19 19:41:14

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

   
    vec2 uv = fragCoord/RENDERSIZE.xy;
    
   
    


       // GAUSSIAN BLUR

    vec2 radius = vec2(6.0);   // PHOTOSHOP HIGHPASS RADIUS +1   vec(2.0)=photoshop's radius 1;
    
    float max = sqrt(radius.x * radius.x + radius.y * radius.y);

    vec3 blur = vec3(0.0);

    
    float sum = 0.0;

    for(float u = -radius.x; u<=radius.x; u++){
        for(float v = -radius.y; v<=radius.y; v++){
           
            float weight = max - sqrt(u * u + v * v)*2;
           
            blur += weight * texture( syn_UserImage, uv + (vec2(u, v)/RENDERSIZE.xy)).xyz;
           
            sum += weight;
        }
    }
   
    blur /= sum;
    

        // 

      
     
        vec3 col=blur;
   		
        vec3 col_orig=(texture(syn_UserImage,uv).rgb);
        
        fragColor = vec4(vec3((col_orig-col))+0.5,1.);
        
       
   
    
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}