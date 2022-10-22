



vec3 rotate(vec3 p, float angle, vec3 axis) {
	vec3 a = normalize(axis);
	float s = sin(angle);
	float c = cos(angle);
	float r = 1.0 - c;
	mat3 m = mat3(
		a.x * a.x * r + c,
		a.y * a.x * r + a.z * s,
		a.z * a.x * r - a.y * s,
		a.x * a.y * r - a.z * s,
		a.y * a.y * r + c,
		a.z * a.y * r + a.x * s,
		a.x * a.z * r + a.y * s,
		a.y * a.z * r - a.x * s,
		a.z * a.z * r + c
	);
	return m * p;
}

float random (in float x) {
	return fract(sin(x) * 1e4);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec3 p) {
	const vec3 step = vec3(110.0, 241.0, 171.0);

	vec3 i = floor(p);
	vec3 f = fract(p);

    // For performance, compute the base input to a
    // 1D random from the integer part of the
    // argument and the incremental change to the
    // 1D based on the 3D -> 1D wrapping
	float n = dot(i, step);

	vec3 u = f * f * (3.0 - 2.0 * f);

	return mix( mix(mix(random(n + dot(step, vec3(0,0,0))),
                        random(n + dot(step, vec3(1,0,0))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,0))),
                        random(n + dot(step, vec3(1,1,0))),
                        u.x),
                u.y),
                mix(mix(random(n + dot(step, vec3(0,0,1))),
                        random(n + dot(step, vec3(1,0,1))),
                        u.x),
                    mix(random(n + dot(step, vec3(0,1,1))),
                        random(n + dot(step, vec3(1,1,1))),
                        u.x),
                u.y),
                u.z);
}

float marble(vec3 p) {
	float m = noise(vec3(10.0 * p.x, 30.0 * (p.x + 0.25 * p.y + p.z), 10.0 * p.z));
	return m * m;
}

vec3 lerp(vec3 a, vec3 b, float f) {
	return a + f * (b - a);
}

float plintoDist(vec3 p) {
	float l1 = max(length(max(abs(mod(p.xz, 4.0) - 2.0) - 0.35, 0.0)), p.y + 1.95);
	float l2 = max(length(max(abs(mod(p.xz, 4.0) - 2.0) - 0.30, 0.0)), p.y + 1.40);
	float l3 = length(vec3(mod(p.xz, 4.0) - 2.0, 4.0 * (p.y + 1.35))) - 0.25;
	return min(min(l1, l2), l3);
}

float colonnaDist(vec3 p) {
	return length(mod(p.xz, 4.0) - 2.0) - 0.2;
}

float capitelloDist(vec3 p) {
	float l1 = length(vec3(mod(p.xz, 4.0) - 2.0, 4.0 * (2.0 - p.y))) - 0.35;
	float l2 = length(vec3(mod(p.xz, 4.0) - 2.0, 4.0 * (1.7 - p.y))) - 0.3;
	float l3 = max(length(max(abs(mod(p.xz, 4.0) - 2.0) -0.21, 0.0)), 1.7 - p.y);
	return min(min(l1, l2), l3);
}

float soffittoDist(vec3 p) {
	float voltax = 1.8 - length(vec2(mod(p.z - 2.0, 4.0) - 2.0, p.y - 2.0));
	float voltaz = 1.8 - length(vec2(mod(p.x - 2.0, 4.0) - 2.0, p.y - 2.0));
	float ceiling = 2.0 - p.y;
	return max(ceiling, max(voltax, voltaz));
}

vec3 lampadaVect(vec3 p) {
	vec3 v = vec3(-(mod(p.xz + 6., 12.0) - 6.0), 1.5 - p.y);
	return vec3(v.x, v.z, v.y);
}

float lampadaGamboDist(vec3 p) {
	return max(length(mod(p.xz - 6.0, 12.0) - 6.0) - 0.01, 1.7 - p.y);
}

float lampadaGloboDist(vec3 p) {
	return length(lampadaVect(p)) - 0.2;
}

float lampadaDist(vec3 p) {
	return min(lampadaGloboDist(p), lampadaGamboDist(p));
}
    
// distance estimator
float sceneDist(vec3 p) {
	return min(min(plintoDist(p), min(colonnaDist(p), soffittoDist(p))), min(capitelloDist(p), lampadaDist(p)));
}

float illuminazione(vec3 rp, vec3 norm) {
	vec3 lv = lampadaVect(rp);
	float l1 = length(lv) - 0.19; // tolta la dimensione della lampada
	float v = max(dot(normalize(lv), norm), 0.1); // 0.1 = luce ambiente
	return v * 3.0 / l1; // intensità della luce inversamente proporzionale alla distanza
}

    
// --------------------------------------------------------------------------------------------------------------------

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = (fragCoord * 2. - RENDERSIZE.xy) / min(RENDERSIZE.x, RENDERSIZE.y);
	vec2 mo = _mouse.xy / RENDERSIZE.xy * 2.0 - 1.0;

	vec3 cameraPos = vec3(0., 0., 3. * TIME);
    
	// angolo visuale
	vec2 a = 0.8 * uv;
	vec3 rayDirection = normalize(vec3(sin(a.x) * cos(a.y), sin(a.y), cos(a.x) * cos(a.y)));
    
	// rotazione vista
	rayDirection = rotate(rayDirection, 0.25 * PI * mo.y, vec3(1., 0., 0.));
	rayDirection = rotate(rayDirection, PI * mo.x, vec3(0., 1., 0.));
    rayDirection.yz = _rotate(rayDirection.yz, lookXY.y*PI+0.75);
    rayDirection.xy = _rotate(rayDirection.xy, lookXY.x*PI);
    rayDirection.xz = _rotate(rayDirection.xz, lookZ*PI);

	const vec3 bgcol = vec3(0.5); // colore foschia
	const vec3 lightcol = vec3(0.95, 0.9, 0.90); // colore luci
    
	vec3 col = bgcol;
	float depth = 0.0;
	float depthPav = 1000.0; // distanza dal pavimento
	bool mirror = false; // riflesso sul pavimento
	vec3 colPav = vec3(1.);

	const float threshold = 0.0001; // soglia distanza minima

	// ray marching
	for (int i = 0; i < 400; i++) {
		vec3 rayPos = cameraPos + rayDirection * depth;
		vec3 rayPosPav;
        
		// riflessione sul pavimento
		if(rayPos.y < -2.) {
			rayPos.y = -4. - rayPos.y;
			if (!mirror) {
				mirror = true;
				depthPav = (cameraPos.y + 2.) / -rayDirection.y;
				rayPosPav = cameraPos + rayDirection * depthPav;
				int b1 = int(mod(rayPosPav.x + rayPosPav.z, 2.0));
				int b2 = int(mod(rayPosPav.x - rayPosPav.z, 2.0));
				int b3 = (b1 > 0) ? b2 : 1 - b2;
				// colore pavimento
				colPav = vec3(0.20 + 0.78 * float(b3), 0.22 + 0.78 * float(b3), 0.20 + 0.76 * float(b3));
			}
		}
        
		float dist = sceneDist(rayPos);
		if (dist < threshold) {
            
            // soffitto
			// vettore normale alla superficie (molto approssimato)
			vec3 norm = vec3(0.0, -1.0, 0.0);
			col = vec3(1.0, 0.9, 0.8);

			// plinto
			if(plintoDist(rayPos)  < threshold) {
				col = lerp(vec3(0.45, 0.50, 0.55), vec3(0.25, 0.30, 0.35), marble(rayPos));
				vec2 v = mod(rayPos.xz, 4.0) - 2.0;
				norm = normalize(vec3(v.x, 0.0, v.y));
			}
			// colonna
			else if(colonnaDist(rayPos) < threshold) {
				col = lerp(vec3(0.90, 0.90, 0.85), vec3(0.80, 0.60, 0.55), marble(rayPos));
				vec2 v = mod(rayPos.xz, 4.0) - 2.0;
				norm = normalize(vec3(v.x, 0.0, v.y));
			}
			// capitello
			else if(capitelloDist(rayPos) < threshold) {
				col = lerp(vec3(0.45, 0.50, 0.55), vec3(0.25, 0.30, 0.35), marble(rayPos));
				vec2 v = mod(rayPos.xz, 4.0) - 2.0;
				norm = normalize(vec3(v.x, 0.0, v.y));
			}
			// lampada
			else if (lampadaGamboDist(rayPos) < threshold) col = vec3(0.0);
			else if (lampadaGloboDist(rayPos) < threshold) col = lightcol+highhits*flashy;

			// luci
			col = col * (lightcol+highhits*flashy) * illuminazione(rayPos, norm);
			if (mirror) colPav = colPav * (lightcol+highhits*flashy) * illuminazione(rayPosPav, vec3(0.0, 0.3, 0.0));
            
			col /= 0.6 + 0.03 * float(i); // accentua i bordi
			break;
        }
		depth += dist;
    }

	// colore riflesso nel pavimento
	if (mirror) col = 0.5 * col + 0.5 * colPav;

	// fog
	col = lerp(col, bgcol, clamp(1.0 - 10.0 / (min(depthPav, depth) + 5.0), 0., 1.));
    
    // Output to screen
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}