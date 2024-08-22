
convertStatePartials(state::AbstractArray, ::Type{MEE}, ::Type{Cartesian}, mu) = cartWrtMee(state, mu)
function cartWrtMee(mee::AbstractArray, mu)
    p=mee[1];f=mee[2];g=mee[3];h=mee[4];k=mee[5];L=mee[6]
    t2 = cos(L);
    t3 = sin(L);
    t4 = f*h;
    t5 = g*k;
    t6 = h*h;
    t7 = h*h*h;
    t8 = k*k;
    t9 = k*k*k;
    t18 = 1.0/p;
    t26 = sqrt(mu);
    t10 = f*t2;
    t11 = f+t2;
    t12 = h*t2;
    t13 = k*t2;
    t14 = g*t3;
    t15 = g+t3;
    t16 = h*t3;
    t17 = k*t3;
    t20 = h*t4;
    t21 = g*t6;
    t22 = f*t8;
    t23 = k*t5;
    t24 = k*t4*2.0;
    t25 = h*t5*2.0;
    t27 = -t6;
    t28 = -t8;
    t32 = t2*t6;
    t36 = t3*t8;
    t42 = sqrt(t18);
    t47 = t6+t8+1.0;
    t29 = -t16;
    t34 = k*t12*2.0;
    t37 = k*t16*2.0;
    t44 = t2*t28;
    t45 = t3*t27;
    t48 = 1.0/t42;
    t49 = t10+t14+1.0;
    t50 = t8+t27+1.0;
    t51 = t6+t28+1.0;
    t52 = 1.0/t47;
    t57 = t4+t5+t12+t17;
    t53 = t52*t52;
    t54 = t13+t29;
    t55 = 1.0/t49;
    t58 = t2+t32+t37+t44;
    t59 = t3+t34+t36+t45;
    t60 = h*k*t26*t42*t52*2.0;
    t56 = t55*t55;
    J00 = t52*t55*t58;
    J01 = -p*t2*t52*t56*t58;
    J02 = -p*t3*t52*t56*t58;
    J03 = k*p*t53*t55*t59*2.0;
    J04 = p*t53*t55*(t13+t54-t3*t7+t6*t13*2.0+t8*t16)*-2.0;
    J05 = p*t52*t56*(-t15-t21+t23+t24+t34+t36+t45);
    J10 = t52*t55*t59;
    J11 = -p*t2*t52*t56*t59;
    J12 = -p*t3*t52*t56*t59;
    J13 = p*t53*t55*(-t13+t16*2.0-t2*t9+t6*t13+t8*t16*2.0)*-2.0;
    J14 = h*p*t53*t55*t58*2.0;
    J15 = -p*t52*t56*(-t11+t20-t22+t25+t32+t37+t44);
    J20 = t52*t54*t55*-2.0;
    J21 = p*t2*t52*t54*t56*2.0;
    J22 = p*t3*t52*t54*t56*2.0;
    J23 = p*t53*t55*t59*2.0;
    J24 = p*t53*t55*t58*-2.0;
    J25 = p*t52*t56*t57*2.0;
    J30 = (t15*t26*1.0/(t48*t48*t48)*t51*t52)/2.0-h*k*t11*t26*1.0/(t48*t48*t48)*t52;
    J31 = t60;
    J32 = -t26*t42*t51*t52;
    J33 = k*t26*t42*t53*(-t11+t20-t22+t25+t32+t37+t44)*-2.0;
    J34 = t26*t42*t53*(t5+t17+t57+t2*t7+t4*t6+t5*t6*2.0+t6*t17*2.0+t4*t28+t12*t28)*2.0;
    J35 = -t26*t42*t52*t58;
    J40 = t11*t26*1.0/(t48*t48*t48)*t50*t52*(-1.0/2.0)+h*k*t15*t26*1.0/(t48*t48*t48)*t52;
    J41 = t26*t42*t50*t52;
    J42 = -t60;
    J43 = t26*t42*t53*(t4+t12+t57+t3*t9+t4*t8*2.0+t5*t8+t8*t12*2.0+t5*t27+t17*t27)*-2.0;
    J44 = h*t26*t42*t53*(-t15-t21+t23+t24+t34+t36+t45)*2.0;
    J45 = -t26*t42*t52*t59;
    J50 = -t26*1.0/(t48*t48*t48)*t52*t57;
    J51 = h*t26*t42*t52*2.0;
    J52 = k*t26*t42*t52*2.0;
    J53 = t26*t42*t53*(-t11+t20-t22+t25+t32+t37+t44)*-2.0;
    J54 = t26*t42*t53*(-t15-t21+t23+t24+t34+t36+t45)*-2.0;
    J55 = t26*t42*t52*t54*2.0;

    return SA[
        J00 J01 J02 J03 J04 J05;
        J10 J11 J12 J13 J14 J15;
        J20 J21 J22 J23 J24 J25;
        J30 J31 J32 J33 J34 J35;
        J40 J41 J42 J43 J44 J45;
        J50 J51 J52 J53 J54 J55
    ]
end

function meeWrtCart(cart::AbstractArray, mu)
    # Put variables in correct form
    rx = cart[1]
    ry = cart[2]
    rz = cart[3]
    vx = cart[4]
    vy = cart[5]
    vz = cart[6]

    t2 = rx*ry;
    t3 = rx*vy;
    t4 = ry*vx;
    t5 = rx*vz;
    t6 = rz*vx;
    t7 = ry*vz;
    t8 = rz*vy;
    t9 = rz*vz;
    t10 = rx*rx;
    t13 = ry*ry;
    t16 = rz*rz;
    t19 = vx*vx;
    t20 = vx*vx*vx;
    t22 = vy*vy;
    t23 = vy*vy*vy;
    t25 = vz*vz;
    t26 = vz*vz*vz;
    t32 = 1.0/mu;
    t27 = t25*t25;
    t28 = t2*vx;
    t29 = t2*vy;
    t30 = rz*t5;
    t31 = rz*t7;
    t33 = -t4;
    t34 = -t6;
    t35 = -t8;
    t36 = rx*t3;
    t37 = ry*t4;
    t38 = rz*t6;
    t39 = rz*t8;
    t42 = rx*t6*2.0;
    t44 = ry*t8*2.0;
    t46 = t3*vx*2.0;
    t47 = t5*vx*2.0;
    t48 = t4*vy*2.0;
    t49 = t7*vy*2.0;
    t50 = t6*vz*2.0;
    t51 = t8*vz*2.0;
    t53 = t5*t6*2.0;
    t54 = t7*t8*2.0;
    t55 = t3*vy*2.0;
    t56 = t4*vx*2.0;
    t59 = t5*vz*2.0;
    t60 = t6*vx*2.0;
    t61 = rx*t5*2.0;
    t63 = t7*vz*2.0;
    t64 = t8*vy*2.0;
    t65 = ry*t7*2.0;
    t77 = t3*t3;
    t78 = t4*t4;
    t79 = t5*t5;
    t80 = t6*t6;
    t81 = t7*t7;
    t82 = t8*t8;
    t86 = mu*t2*t9*2.0;
    t102 = mu*t2*t4*2.0;
    t103 = mu*t2*t3*2.0;
    t104 = mu*t2*t25*2.0;
    t105 = mu*t2*t6*vz*-2.0;
    t109 = mu*t2*t7*vz;
    t112 = mu*t2*t5*vz;
    t125 = t10+t13+t16;
    t127 = mu*ry*t7*t8*-2.0;
    t40 = t28*2.0;
    t41 = t29*2.0;
    t43 = t30*2.0;
    t45 = t31*2.0;
    t57 = t36*2.0;
    t58 = t37*2.0;
    t62 = t38*2.0;
    t66 = t39*2.0;
    t67 = -t28;
    t71 = -t31;
    t74 = -t48;
    t76 = -t50;
    t83 = t28*vy*-2.0;
    t84 = -t53;
    t85 = -t54;
    t87 = mu*t28*vy*4.0;
    t88 = ry*t33;
    t89 = -t56;
    t91 = rz*t34;
    t93 = -t63;
    t95 = mu*t2*t51;
    t96 = mu*t3*t36;
    t97 = mu*t4*t37;
    t99 = mu*rx*t79;
    t100 = mu*t82;
    t101 = mu*ry*t81;
    t108 = mu*t3*t29;
    t111 = mu*t4*t38;
    t114 = t3+t33;
    t115 = t5+t34;
    t116 = t7+t35;
    t119 = mu*t5*t42;
    t122 = mu*t6*t34;
    t124 = mu*t8*t35;
    t126 = mu*t4*t29*-2.0;
    t131 = mu*t28*t33;
    t133 = -t109;
    t134 = mu*rz*t3*t35;
    t147 = 1.0/t125;
    t150 = sqrt(t125);
    t52 = t40*vy;
    t68 = -t40;
    t72 = -t45;
    t90 = -t58;
    t92 = -t62;
    t113 = ry*t100;
    t118 = mu*t3*t40;
    t121 = -t96;
    t123 = -t99;
    t128 = t114*t114;
    t129 = t115*t115;
    t130 = t116*t116;
    t132 = rx*t122;
    t135 = rx*t114*2.0;
    t136 = ry*t114*2.0;
    t137 = rx*t115*2.0;
    t138 = rz*t115*2.0;
    t139 = ry*t116*2.0;
    t140 = rz*t116*2.0;
    t141 = t114*vx*2.0;
    t142 = t114*vy*2.0;
    t143 = t115*vx*2.0;
    t144 = t115*vz*2.0;
    t145 = t116*vy*2.0;
    t146 = t116*vz*2.0;
    t151 = 1.0/t150;
    t153 = (t3*t3*t3)*t150*vy;
    t155 = (t5*t5*t5)*t150*vz;
    t156 = (t7*t7*t7)*t150*vz;
    t157 = t2*t27*t150*2.0;
    t159 = t2*t7*t26*t150;
    t160 = t3*t82*t150*vy;
    t162 = t2*t5*t26*t150;
    t163 = (t3*t3*t3)*t150*vx;
    t165 = t6*t8*t22*t150;
    t166 = t80*t150*vx*vy;
    t168 = (t8*t8*t8)*t150*vz;
    t169 = t2*t9*t25*t150*2.0;
    t170 = t3*t9*t25*t150*2.0;
    t171 = t3*t51*t150*vy;
    t172 = t4*t9*t25*t150*2.0;
    t173 = t4*t50*t150*vx;
    t174 = t6*t8*t9*t150*2.0;
    t175 = t3*t6*t8*t150*vy;
    t176 = t3*t80*t150*vx;
    t179 = t3*t5*t9*t150*vz;
    t180 = t8*t77*t150*vz;
    t183 = t3*t79*t150*vx;
    t186 = t8*t80*t150*vz;
    t187 = t2*t6*t26*t150*2.0;
    t188 = t2*t19*t50*t150;
    t189 = t2*t8*t26*t150*2.0;
    t191 = t25*t28*t150*vy*6.0;
    t192 = t3*t6*t150*vx*vz*4.0;
    t193 = t4*t8*t150*vy*vz*4.0;
    t194 = t33*t78*t150*vx;
    t196 = t6*t29*t150*vz*8.0;
    t199 = t77*t150*vx*vy*3.0;
    t200 = t78*t150*vx*vy*3.0;
    t203 = t22*t80*t150;
    t204 = t4*t19*t29*t150*3.0;
    t205 = t3*t22*t28*t150*3.0;
    t206 = t2*t19*t22*t150*6.0;
    t207 = t2*t19*t25*t150*2.0;
    t208 = t2*t22*t25*t150*2.0;
    t209 = t6*t79*t150*vz*3.0;
    t210 = t3*(t9*t9)*t150*2.0;
    t211 = t4*(t9*t9)*t150*2.0;
    t212 = t9*t77*t150*3.0;
    t213 = t9*t78*t150*3.0;
    t214 = t8*t81*t150*vz*3.0;
    t215 = t3*t5*t150*vx*vz*3.0;
    t216 = t4*t7*t150*vy*vz*3.0;
    t217 = t6*t8*t25*t150*2.0;
    t219 = t2*t8*t22*t150*vz*-2.0;
    t224 = t2*t6*t22*t150*vz*4.0;
    t225 = t6*t28*t150*vy*vz*4.0;
    t226 = t33*t80*t150*vx;
    t228 = t33*t78*t150*vy;
    t229 = t34*t80*t150*vz;
    t231 = t59*t77*t150;
    t237 = t33*t82*t150*vy;
    t238 = t33*t80*t150*vy;
    t240 = t7*t9*t33*t150*vz;
    t241 = t4*t6*t33*t150*vz;
    t242 = t33*t81*t150*vy;
    t244 = t34*t82*t150*vz;
    t245 = t4*t25*t28*t150;
    t246 = t3*t80*t150*vy;
    t247 = t2*t5*t19*t150*vz;
    t248 = t2*t7*t22*t150*vz;
    t250 = t3*t25*t29*t150;
    t251 = t3*t8*t9*t150*vz;
    t255 = t4*t25*t29*t150*3.0;
    t256 = t3*t25*t28*t150*3.0;
    t259 = t6*t77*t150*vz*3.0;
    t260 = t4*t9*t51*t150;
    t261 = t7*t8*t48*t150;
    t262 = t8*t78*t150*vz*3.0;
    t270 = t7*t78*t150*vz*-2.0;
    t271 = t4*t22*t28*t150*3.0;
    t272 = t3*t19*t29*t150*3.0;
    t273 = t5*t80*t150*vz*3.0;
    t274 = t7*t82*t150*vz*3.0;
    t276 = t3*t6*t9*t150*vz*-2.0;
    t277 = t3*t5*t6*t150*vx*-2.0;
    t284 = t6*t8*t33*t150*vy;
    t286 = t6*t9*t33*t150*vz;
    t148 = -t140;
    t149 = -t146;
    t152 = t151*t151*t151;
    t158 = -t157;
    t195 = -t156;
    t197 = t2*t25*t81*t151;
    t198 = t2*t25*t79*t151;
    t201 = -t169;
    t202 = -t174;
    t218 = -t189;
    t220 = -t191;
    t221 = -t192;
    t222 = -t193;
    t223 = -t196;
    t227 = -t162;
    t230 = (t2*t2)*t27*t151;
    t233 = -t205;
    t234 = -t206;
    t235 = -t207;
    t236 = -t208;
    t239 = -t209;
    t243 = -t217;
    t263 = -t224;
    t266 = t136+t138;
    t267 = t137+t139;
    t268 = t142+t144;
    t269 = t143+t145;
    t275 = -t256;
    t278 = -t259;
    t281 = t2*t6*t8*t25*t151*2.0;
    t282 = -t272;
    t283 = -t247;
    t285 = -t250;
    t287 = -t274;
    t288 = t2*t8*t22*t34*t151;
    t289 = t6*t28*t34*t151*vy;
    t291 = (t2*t2)*t19*t25*t151;
    t292 = (t2*t2)*t22*t25*t151;
    t293 = t25*t28*t29*t151*3.0;
    t294 = (t2*t2)*t19*t22*t151*3.0;
    t295 = t2*t6*t22*t34*t151;
    t296 = t128+t129+t130;
    t264 = -t197;
    t265 = -t198;
    t279 = t135+t148;
    t280 = t141+t149;
    t297 = 1.0/t296;
    t298 = sqrt(t296);
    t299 = 1.0/t298;
    t300 = rx*t298;
    t301 = ry*t298;
    t302 = t298*vx;
    t303 = t298*vy;
    t304 = t3*t298;
    t305 = t4*t298;
    t306 = mu*t2*t298;
    t307 = mu*t9*t298;
    t308 = mu*t16*t298;
    t309 = mu*t28*t298;
    t310 = mu*t29*t298;
    t311 = mu*t30*t298;
    t312 = mu*t31*t298;
    t313 = t33*t298;
    t314 = mu*t36*t298;
    t319 = mu*t88*t298;
    t320 = mu*t91*t298;
    t321 = mu*rz*t35*t298;
    t322 = t114+t298;
    t332 = t6*t8*t150*t298*2.0;
    t337 = t79*t150*t298;
    t339 = t81*t150*t298;
    t356 = t22*t28*t150*t298*-2.0;
    t357 = t19*t29*t150*t298*-2.0;
    t360 = t6*t41*t151*t298*vz;
    t318 = -t314;
    t323 = 1.0/t322;
    t328 = t28*t150*t303*4.0;
    t329 = t50*t150*t303;
    t330 = t46*t150*t303;
    t335 = t77*t150*t303;
    t336 = t78*t150*t302;
    t338 = t80*t150*t302;
    t340 = t82*t150*t303;
    t341 = t4*t150*t302*vy*-2.0;
    t342 = t77*t150*t302;
    t343 = t79*t150*t302;
    t344 = t78*t150*t303;
    t345 = t5*t150*t304*vz;
    t346 = t7*t150*t305*vz;
    t347 = t6*t8*t150*t303;
    t348 = t80*t150*t303;
    t349 = t81*t150*t303;
    t354 = t6*t150*t304*vz*-2.0;
    t355 = t8*t150*t305*vz*-2.0;
    t358 = t5*t6*t150*t302*-2.0;
    t359 = t7*t8*t150*t303*-2.0;
    t361 = t36+t39+t67+t71+t300;
    t362 = t29+t30+t88+t91+t301;
    t367 = (t266*t299)/2.0;
    t368 = (t268*t299)/2.0;
    t369 = (t279*t299)/2.0;
    t370 = (t280*t299)/2.0;
    t375 = (t3*t267*t299)/2.0;
    t376 = (t4*t267*t299)/2.0;
    t379 = (t3*t269*t299)/2.0;
    t380 = (t4*t269*t299)/2.0;
    t383 = t4*t266*t299*(-1.0/2.0);
    t387 = t4*t268*t299*(-1.0/2.0);
    t391 = t4*t279*t299*(-1.0/2.0);
    t392 = t4*t280*t299*(-1.0/2.0);
    t395 = t77+t78+t79+t80+t81+t82+t83+t84+t85+t304+t313;
    t324 = t323*t323;
    t325 = rz*t323;
    t326 = t323*vz;
    t363 = t361*t361;
    t364 = t362*t362;
    t365 = 1.0/t361;
    t371 = ry+t367;
    t372 = t368+vy;
    t373 = t3*t367;
    t374 = t4*t367;
    t377 = t3*t368;
    t378 = t4*t368;
    t381 = rx+t369;
    t382 = t370+vx;
    t384 = t3*t369;
    t385 = t4*t369;
    t388 = t3*t370;
    t389 = t4*t370;
    t396 = 1.0/t395;
    t406 = t97+t101+t105+t108+t111+t112+t113+t126+t127+t163+t168+t175+t176+t179+t180+t183+t186+t187+t188+t194+t195+t204+t214+t226+t227+t255+t262+t263+t270+t276+t277+t282+t283+t284+t285+t286+t287+t310+t311+t319+t320+t336+t338+t342+t343+t346+t347+t355+t357+t358;
    t407 = t95+t118+t119+t121+t123+t131+t132+t133+t134+t153+t155+t159+t160+t218+t219+t225+t228+t229+t231+t233+t237+t238+t239+t240+t241+t242+t244+t245+t246+t248+t251+t260+t261+t271+t273+t275+t278+t309+t312+t318+t321+t335+t340+t344+t345+t348+t349+t354+t356+t359;
    t327 = -t325;
    t366 = 1.0/t363;
    t393 = t147*t297*t363;
    t394 = t147*t297*t364;
    t397 = t396*t396;
    t400 = t41+t43+t90+t92+t301+t373+t383;
    t401 = t55+t59+t74+t76+t303+t377+t387;
    t402 = t57+t66+t68+t72+t300+t384+t391;
    t403 = t46+t51+t89+t93+t302+t388+t392;
    t404 = t393+t394;
    t405 = 1.0/t404;
    J00 = t32*t268;
    J01 = -t32*t280;
    J02 = -t32*t269;
    J03 = -t32*t266;
    J04 = t32*t279;
    J05 = t32*t267;
    J10 = t32*t151*t396*(t87+t122+t124+t203+t220+t230+t281+t288+t289+t291+t292+t294+mu*t54-mu*t77*3.0-mu*t79*3.0-mu*t81-mu*t304*2.0+mu*t305+(t77*t77)*t151+(t79*t79)*t151-(t5*t5*t5)*t6*t151*3.0+(t9*t9)*t77*t151+(t3*t3)*t151*t304+mu*t5*t6*4.0+mu*t4*t33+mu*t28*t368+mu*t31*t368-t23*t28*t150*6.0+t22*t77*t150*3.0+t22*t78*t150*3.0+t25*t77*t150*6.0+t22*t81*t150+t25*t78*t150+t22*t82*t150+t25*t79*t150*3.0+t25*t80*t150*3.0+t25*t81*t150+t25*t82*t150+t77*t79*t151*2.0+t77*t80*t151+t77*t82*t151+t79*t80*t151*3.0+t22*t150*t304*2.0-t22*t150*t305*2.0+t25*t150*t304*2.0+t79*t151*t304+t80*t151*t304+t82*t151*t304+t83*t151*t304-(mu*t36*t268*t299)/2.0-(mu*t39*t268*t299)/2.0-t5*t6*t25*t150*6.0-t7*t8*t22*t150*2.0-t7*t8*t25*t150*2.0-t5*t6*t77*t151*3.0+t5*t34*t80*t151-t5*t6*t151*t304*2.0+t4*t28*t151*t303-t28*t77*t151*vy*3.0+t77*t150*t368*vy+t78*t150*t368*vy+t80*t150*t368*vy+t81*t150*t368*vy+t82*t150*t368*vy-t6*t150*t303*vz*2.0+t5*t150*t377*vz-t2*t3*t9*t25*t151*2.0+t3*t8*t9*t34*t151+t2*t9*t25*t33*t151+t22*t67*t150*t268*t299+t4*t28*t33*t151*vy-t3*t5*t28*t151*vz*3.0+t3*t6*t28*t151*vz*4.0-t3*t8*t29*t151*vz*2.0+t4*t8*t41*t151*vz+t6*t28*t33*t151*vz+t7*t29*t33*t151*vz+t7*t29*t151*t298*vz-t8*t29*t151*t298*vz*2.0-t3*t6*t150*vy*vz*6.0+t4*t6*t150*vy*vz*4.0+t7*t35*t150*t268*t299*vy+t3*t34*t150*t268*t299*vz)-rx*t32*t152*t396*t407-t32*t151*t397*t401*t407;
    J11 = -t32*t151*t396*(t104+t158+t165+t166+t170+t171+t172+t173+t199+t200+t215+t216+t221+t222+t234+t235+t236+t243+t264+t265+t293+t295-t307+t330+t341+t360+(t4*t4*t4)*t9*t151-mu*t3*t9*2.0+mu*t2*t19*2.0+mu*t28*t370+mu*t31*t370-mu*t36*vx*2.0-mu*t300*vx+t51*t150*t303-(mu*t36*t280*t299)/2.0-(mu*t39*t280*t299)/2.0+t2*t23*t28*t151*3.0+t2*t22*t54*t151+t2*t25*t54*t151+t6*t8*t78*t151+t4*t9*t80*t151+t4*t9*t81*t151+t4*t9*t82*t151-t2*t22*t77*t151-t2*t22*t78*t151*3.0-t2*t25*t77*t151*2.0-t2*t22*t81*t151-t2*t25*t80*t151*3.0+t6*t8*t151*t313-t2*t22*t151*t304+t2*t22*t151*t305*2.0-t2*t25*t151*t304+t7*t44*t151*t303+t33*t37*t151*t303+t37*t78*t151*vy+t37*t81*t151*vy+t37*t82*t151*vy-t81*t151*t301*vy+t77*t150*t370*vy+t78*t150*t370*vy+t80*t150*t370*vy+t81*t150*t370*vy+t82*t150*t370*vy-t7*t150*t303*vz*2.0+t5*t150*t388*vz-t4*t7*t8*t9*t151*2.0+t2*t5*t6*t25*t151*3.0+t2*t4*t25*t33*t151+t2*t8*t22*t35*t151+t2*t8*t25*t35*t151+t22*t67*t150*t280*t299-t7*t8*t37*t151*vy*2.0+t8*t35*t151*t301*vy+t3*t6*t29*t151*vz*3.0-t4*t6*t29*t151*vz*4.0+t7*t35*t150*t280*t299*vy+t3*t34*t150*t280*t299*vz)-ry*t32*t152*t396*t407+t32*t151*t397*t403*t407;
    J12 = -t32*t151*t396*(t4*(t8*t8*t8)*t151+(t4*t4*t4)*t8*t151-(t5*t5*t5)*t9*t151+(t6*t6*t6)*t9*t151+(t3*t3*t3)*t35*t151+mu*t3*t8*2.0-mu*t7*t298+mu*t8*t298*2.0+mu*t42*vx-mu*t29*vz*2.0+t26*t41*t150+(mu*t28*t269*t299)/2.0+(mu*t31*t269*t299)/2.0-(mu*t36*t269*t299)/2.0-(mu*t39*t269*t299)/2.0-mu*rx*t5*vx*2.0-t3*t8*t22*t150*2.0-t4*t7*t22*t150*2.0+t4*t8*t22*t150*2.0-t3*t8*t25*t150*2.0+t4*t7*t25*t150-t4*t8*t25*t150*4.0-t5*t9*t77*t151*2.0+t4*t8*t80*t151+t6*t9*t77*t151*3.0-t4*t7*t82*t151*2.0+t4*t8*t81*t151+t6*t9*t78*t151-t5*t9*t80*t151*3.0+t6*t9*t79*t151*3.0+t6*t9*t82*t151+t3*t35*t82*t151-t5*t9*t151*t304+t6*t9*t151*t304*2.0+t7*t22*t150*t298*2.0-t8*t22*t150*t298*2.0+t8*t33*t151*t305+t8*t54*t151*t298+t35*t77*t151*t298+t35*t81*t151*t298+t35*t82*t151*t298+t6*t48*t150*vx-t6*t150*t302*vy*2.0+t2*t23*t150*vz*2.0-t19*t29*t150*vz*4.0-t29*t80*t151*vz*4.0+t41*t82*t151*vz+t46*t150*t298*vz+t5*t150*t379*vz+t77*t150*vx*vz*3.0+t78*t150*vx*vz+t79*t150*vx*vz*3.0+t80*t150*vx*vz*3.0+t4*t7*(t9*t9)*t151-t4*t8*(t9*t9)*t151*2.0+t3*(t9*t9)*t35*t151+t2*t3*t6*t22*t151*3.0-t2*t4*t6*t22*t151*3.0+t2*t3*t6*t25*t151*3.0-t2*t7*t9*t25*t151+t2*t8*t9*t25*t151*2.0+t3*t6*t8*t34*t151+t2*t6*t25*t33*t151+t2*t6*t22*t151*t298*2.0+t6*t8*t34*t151*t298+t22*t67*t150*t269*t299+(t77*t150*t269*t299*vy)/2.0+(t78*t150*t269*t299*vy)/2.0+(t80*t150*t269*t299*vy)/2.0+(t81*t150*t269*t299*vy)/2.0+(t82*t150*t269*t299*vy)/2.0+t7*t29*t35*t151*vz-t3*t6*t150*vx*vy*2.0-t5*t6*t150*vx*vz*6.0+t6*t8*t150*vy*vz*3.0+t7*t35*t150*t269*t299*vy+t3*t34*t150*t269*t299*vz)+t32*t151*t397*t407*(t47+t49-t60-t64+t379-t380)-rz*t32*t152*t396*t407;
    J13 = -t32*t151*t396*(t102-t103+t212+t213+t223-t306-t332-mu*rx*t30*2.0+mu*rx*t62+mu*t28*t367+mu*t31*t367+t9*t79*t150*3.0+t9*t80*t150*3.0+t9*t81*t150+t9*t82*t150+t9*t150*t304*2.0-t37*t150*t303*2.0-(mu*t36*t266*t299)/2.0-(mu*t39*t266*t299)/2.0-t3*t6*t8*t150*2.0+t4*t6*t8*t150*3.0-t5*t6*t9*t150*6.0-t7*t8*t9*t150*2.0+t2*t3*t22*t150*3.0-t2*t4*t22*t150*6.0+t2*t3*t25*t150*3.0-t2*t4*t25*t150*2.0+t2*t22*t150*t298*2.0+ry*t81*t150*vy+ry*t82*t150*vy+t4*t37*t150*vy*3.0+t77*t150*t367*vy+t78*t150*t367*vy+t80*t150*t367*vy+t81*t150*t367*vy+t82*t150*t367*vy+t5*t150*t373*vz+t22*t67*t150*t266*t299-ry*t7*t8*t150*vy*2.0+t7*t35*t150*t266*t299*vy+t3*t34*t150*t266*t299*vz)+t32*t151*t397*t400*t407;
    J14 = t32*t151*t396*(t86+t201+t202+t210+t211-t308-t328+t337+t339+(t3*t3*t3)*t150*4.0+mu*rx*t40-mu*t3*t10*2.0-mu*t3*t16*2.0-mu*t10*t298+mu*t28*t369+mu*t31*t369+t3*t79*t150*4.0+t3*t80*t150*2.0+t3*t82*t150*4.0-t4*t82*t150*3.0+t33*t78*t150+t33*t80*t150+t33*t81*t150+t77*t150*t298*3.0+t78*t150*t298+t80*t150*t298+t82*t150*t298*3.0-(mu*t36*t279*t299)/2.0-(mu*t39*t279*t299)/2.0-t3*t5*t6*t150*6.0+t4*t7*t8*t150*4.0-t5*t6*t150*t298*2.0-t7*t8*t150*t298*4.0-t3*t28*t150*vy*9.0+t4*t28*t150*vy*6.0+t77*t150*t369*vy+t78*t150*t369*vy+t80*t150*t369*vy+t81*t150*t369*vy+t82*t150*t369*vy-t5*t28*t150*vz*3.0+t6*t28*t150*vz*4.0-t8*t29*t150*vz*6.0+t7*t41*t150*vz+t5*t150*t384*vz+t22*t67*t150*t279*t299+t7*t35*t150*t279*t299*vy+t3*t34*t150*t279*t299*vz)-t32*t151*t397*t402*t407;
    J15 = t32*t151*t396*((t5*t5*t5)*t150*4.0+mu*rz*t301-mu*t2*t7*2.0+mu*t2*t8*2.0-mu*t5*t10*2.0+mu*t6*t10*2.0+t5*t77*t150*4.0-t6*t77*t150*3.0+t5*t80*t150*6.0-t6*t79*t150*9.0+t34*t80*t150+t34*t82*t150+t5*t150*t304*2.0-t6*t150*t304*2.0+t49*t150*t301+(mu*t28*t267*t299)/2.0+(mu*t31*t267*t299)/2.0-(mu*t36*t267*t299)/2.0-(mu*t39*t267*t299)/2.0+t3*t8*t9*t150*2.0-t4*t7*t9*t150*3.0+t4*t8*t9*t150*4.0+t2*t7*t22*t150*2.0-t2*t8*t22*t150*2.0+t2*t7*t25*t150*4.0-t2*t8*t25*t150*6.0+t4*t6*t33*t150+t6*t28*t150*vy*4.0-t7*t37*t150*vy*2.0+t8*t58*t150*vy-t8*t150*t301*vy*2.0-t3*t28*t150*vz*6.0+t4*t40*t150*vz+t5*t150*t375*vz+t22*t67*t150*t267*t299+(t77*t150*t267*t299*vy)/2.0+(t78*t150*t267*t299*vy)/2.0+(t80*t150*t267*t299*vy)/2.0+(t81*t150*t267*t299*vy)/2.0+(t82*t150*t267*t299*vy)/2.0+t7*t35*t150*t267*t299*vy+t3*t34*t150*t267*t299*vz)+t32*t151*t397*t407*(t42+t44-t61-t65-t375+t376);
    J20 = -t32*t151*t396*(t104+t158+t165+t166+t170+t171+t172+t173+t199+t200+t215+t216+t221+t222+t234+t235+t236+t243+t264+t265+t293+t295+t307+t330+t341+(t3*t3*t3)*t9*t151-mu*t4*t9*2.0+mu*t2*t22*2.0+mu*t29*t368+mu*t30*t368-mu*t37*vy*2.0+mu*t301*vy-(mu*t37*t268*t299)/2.0-(mu*t38*t268*t299)/2.0+t2*t20*t29*t151*3.0+t2*t19*t53*t151+t2*t25*t53*t151+t3*t9*t79*t151+t6*t8*t77*t151+t3*t9*t80*t151+t3*t9*t82*t151-t2*t19*t77*t151*3.0-t2*t19*t79*t151-t2*t25*t77*t151-t2*t25*t78*t151*2.0-t2*t25*t82*t151*3.0+t6*t8*t151*t304-t2*t19*t151*t304*2.0+t2*t19*t151*t305+t2*t25*t151*t305+t3*t36*t151*t302+t36*t77*t151*vx+t36*t79*t151*vx+t36*t80*t151*vx+t79*t151*t300*vx+t80*t151*t300*vx+t77*t150*t368*vx+t78*t150*t368*vx+t79*t150*t368*vx+t80*t150*t368*vx-t6*t150*t302*vz*2.0+t47*t150*t298*vz+t7*t150*t378*vz-t3*t5*t6*t9*t151*2.0+t2*t7*t8*t25*t151*3.0+t2*t4*t19*t33*t151+t2*t6*t19*t34*t151+t2*t6*t25*t34*t151-t19*t29*t150*t268*t299-t5*t6*t36*t151*vx*2.0-t5*t6*t151*t300*vx*2.0+t6*t8*t150*t368*vy-t3*t6*t29*t151*vz*4.0+t4*t6*t29*t151*vz*3.0-t6*t29*t151*t298*vz*2.0+t5*t34*t150*t268*t299*vx+t8*t33*t150*t268*t299*vz)+rx*t32*t152*t396*t406+t32*t151*t397*t401*t406;
    J21 = t32*t151*t396*(t87+t122+t124+t203+t220+t230+t281+t288+t289+t291+t292+t294+t329+mu*t53-mu*t77-mu*t78*3.0-mu*t79-mu*t81*3.0-mu*t304+mu*t305*2.0+(t78*t78)*t151+(t81*t81)*t151-(t7*t7*t7)*t8*t151*3.0+(t9*t9)*t78*t151+mu*t7*t8*4.0+mu*t29*t370+mu*t30*t370-t20*t29*t150*6.0+t19*t77*t150*3.0+t19*t78*t150*3.0+t19*t79*t150+t19*t80*t150+t25*t77*t150+t25*t78*t150*6.0+t25*t79*t150+t25*t80*t150+t25*t81*t150*3.0+t25*t82*t150*3.0+t78*t80*t151+t78*t81*t151*2.0+t78*t82*t151+t81*t82*t151*3.0+t19*t150*t304*2.0-t19*t150*t305*2.0-t25*t150*t305*2.0+t52*t151*t305+t54*t151*t305+t78*t151*t313+t80*t151*t313+t81*t151*t313+t82*t151*t313-(mu*t37*t280*t299)/2.0-(mu*t38*t280*t299)/2.0-t5*t6*t19*t150*2.0-t5*t6*t25*t150*2.0-t7*t8*t25*t150*6.0-t7*t8*t78*t151*3.0+t7*t35*t82*t151+t3*t67*t151*t303+t77*t150*t370*vx+t78*t150*t370*vx+t79*t150*t370*vx+t80*t150*t370*vx-t28*t78*t151*vy*3.0+t67*t77*t151*vy+t7*t150*t389*vz-t2*t3*t9*t25*t151-t2*t4*t9*t25*t151*2.0+t6*t8*t9*t33*t151-t19*t29*t150*t280*t299+t6*t8*t150*t370*vy-t4*t6*t28*t151*vz*2.0-t4*t7*t29*t151*vz*3.0+t4*t8*t29*t151*vz*4.0+t3*t6*t40*t151*vz+t3*t29*t35*t151*vz+t3*t5*t67*t151*vz+t6*t40*t151*t298*vz+t5*t67*t151*t298*vz+t3*t6*t150*vy*vz*4.0-t4*t6*t150*vy*vz*6.0+t5*t34*t150*t280*t299*vx+t8*t33*t150*t280*t299*vz)+ry*t32*t152*t396*t406-t32*t151*t397*t403*t406;
    J22 = -t32*t151*t396*(t3*(t6*t6*t6)*t151+(t3*t3*t3)*t6*t151-(t7*t7*t7)*t9*t151+(t8*t8*t8)*t9*t151+(t6*t6*t6)*t33*t151+(t6*t6*t6)*t151*t298+mu*t4*t6*2.0+mu*t5*t298-mu*t6*t298*2.0+mu*t44*vy-mu*t28*vz*2.0+t26*t40*t150-(mu*t29*t269*t299)/2.0-(mu*t30*t269*t299)/2.0+(mu*t37*t269*t299)/2.0+(mu*t38*t269*t299)/2.0-mu*ry*t7*vy*2.0-t3*t5*t19*t150*2.0+t3*t6*t19*t150*2.0-t4*t6*t19*t150*2.0+t3*t6*t22*t150*2.0-t4*t6*t22*t150*2.0+t3*t5*t25*t150-t3*t6*t25*t150*4.0-t4*t6*t25*t150*2.0-t3*t5*t80*t151*2.0+t3*t6*t79*t151+t3*t6*t82*t151-t7*t9*t78*t151*2.0+t8*t9*t77*t151+t8*t9*t78*t151*3.0+t8*t9*t80*t151-t7*t9*t82*t151*3.0+t8*t9*t81*t151*3.0+t6*t33*t78*t151+t6*t33*t82*t151-t5*t19*t150*t298*2.0+t7*t9*t151*t305+t6*t19*t150*t298*2.0-t8*t9*t151*t305*2.0+t6*t22*t150*t298*2.0+t6*t77*t151*t298+t6*t78*t151*t298-t5*t80*t151*t298*2.0+t6*t79*t151*t298+t6*t82*t151*t298+t6*t83*t151*t298+t2*t20*t150*vz*2.0-t22*t28*t150*vz*4.0+t40*t80*t151*vz-t4*t150*t303*vz*2.0+t77*t150*vy*vz+t78*t150*vy*vz*3.0+t80*t150*vy*vz*3.0+t81*t150*vy*vz*3.0+t82*t150*vy*vz*3.0+t3*t5*(t9*t9)*t151-t3*t6*(t9*t9)*t151*2.0+t6*(t9*t9)*t33*t151+t2*t4*t8*t25*t151*3.0-t2*t5*t9*t25*t151+t2*t6*t9*t25*t151*2.0+t2*t3*t25*t35*t151+t19*t29*t150*t269*t299-(t77*t150*t269*t299*vx)/2.0-(t78*t150*t269*t299*vx)/2.0-(t79*t150*t269*t299*vx)/2.0-(t80*t150*t269*t299*vx)/2.0-t3*t6*t28*t151*vy*3.0+t4*t6*t28*t151*vy*3.0-t6*t8*t29*t151*vz*4.0+t5*t28*t34*t151*vz-t7*t8*t150*vy*vz*6.0+t5*t6*t150*t269*t299*vx-(t6*t8*t150*t269*t299*vy)/2.0-(t4*t7*t150*t269*t299*vz)/2.0+t4*t8*t150*t269*t299*vz)-t32*t151*t397*t406*(t47+t49-t60-t64+t379-t380)+rz*t32*t152*t396*t406;
    J23 = t32*t151*t396*(t86+t201+t202+t210+t211+t308+t328-t337-t339-(t3*t3*t3)*t150+(t4*t4*t4)*t150*4.0+mu*ry*t41-mu*t4*t13*2.0-mu*t4*t16*2.0+mu*t13*t298+mu*t29*t367+mu*t30*t367-t3*t79*t150-t3*t80*t150*3.0+t4*t80*t150*4.0+t4*t81*t150*4.0+t4*t82*t150*2.0+t54*t150*t298-t77*t150*t298-t78*t150*t298*3.0-t80*t150*t298*3.0-(mu*t37*t266*t299)/2.0-(mu*t38*t266*t299)/2.0+t3*t5*t6*t150*4.0-t4*t7*t8*t150*6.0+t3*t8*t35*t150+t5*t6*t150*t298*4.0+t8*t35*t150*t298+t77*t150*t367*vx+t78*t150*t367*vx+t79*t150*t367*vx+t80*t150*t367*vx+t3*t28*t150*vy*6.0-t4*t28*t150*vy*9.0-t6*t28*t150*vz*6.0-t7*t29*t150*vz*3.0+t8*t29*t150*vz*4.0+t5*t40*t150*vz+t7*t150*t374*vz-t19*t29*t150*t266*t299+t6*t8*t150*t367*vy+t5*t34*t150*t266*t299*vx+t8*t33*t150*t266*t299*vz)-t32*t151*t397*t400*t406;
    J24 = -t32*t151*t396*(-t102+t103+t212+t213+t223+t306+t332-mu*ry*t31*2.0+mu*ry*t66+mu*t29*t369+mu*t30*t369+t9*t79*t150+t9*t80*t150+t9*t81*t150*3.0+t9*t82*t150*3.0-t9*t150*t305*2.0+t57*t150*t302-(mu*t37*t279*t299)/2.0-(mu*t38*t279*t299)/2.0+t3*t6*t8*t150*3.0-t4*t6*t8*t150*2.0-t5*t6*t9*t150*2.0-t2*t3*t19*t150*6.0-t7*t8*t9*t150*6.0+t2*t4*t19*t150*3.0-t2*t3*t25*t150*2.0+t2*t4*t25*t150*3.0-t2*t19*t150*t298*2.0+rx*t79*t150*vx+rx*t80*t150*vx+t3*t36*t150*vx*3.0+t77*t150*t369*vx+t78*t150*t369*vx+t79*t150*t369*vx+t80*t150*t369*vx+t7*t150*t385*vz-t19*t29*t150*t279*t299-rx*t5*t6*t150*vx*2.0+t6*t8*t150*t369*vy+t5*t34*t150*t279*t299*vx+t8*t33*t150*t279*t299*vz)+t32*t151*t397*t402*t406;
    J25 = -t32*t151*t396*((t7*t7*t7)*t150*-4.0+(t8*t8*t8)*t150+mu*rz*t300+mu*t2*t5*2.0-mu*t2*t6*2.0+mu*t7*t13*2.0-mu*t8*t13*2.0+t36*t47*t150-t7*t78*t150*4.0+t8*t77*t150+t8*t78*t150*3.0+t8*t80*t150-t7*t82*t150*6.0+t8*t81*t150*9.0+t7*t150*t305*2.0-t8*t150*t305*2.0+t47*t150*t300+(mu*t29*t267*t299)/2.0+(mu*t30*t267*t299)/2.0-(mu*t37*t267*t299)/2.0-(mu*t38*t267*t299)/2.0+t3*t5*t9*t150*3.0-t3*t6*t9*t150*4.0-t4*t6*t9*t150*2.0-t2*t5*t19*t150*2.0+t2*t6*t19*t150*2.0-t2*t6*t22*t150*4.0-t2*t5*t25*t150*4.0+t2*t6*t25*t150*6.0-t6*t36*t150*vx*2.0-t6*t150*t300*vx*2.0-t3*t29*t150*vz*2.0+t4*t29*t150*vz*6.0+t7*t150*t376*vz-t19*t29*t150*t267*t299+(t77*t150*t267*t299*vx)/2.0+(t78*t150*t267*t299*vx)/2.0+(t79*t150*t267*t299*vx)/2.0+(t80*t150*t267*t299*vx)/2.0+t5*t34*t150*t267*t299*vx+(t6*t8*t150*t267*t299*vy)/2.0+t8*t33*t150*t267*t299*vz)-t32*t151*t397*t406*(t42+t44-t61-t65-t375+t376);
    J30 = t326-t115*t324*t372;
    J31 = t115*t324*t382;
    J32 = -t323*vx+(t115*t269*t299*t324)/2.0;
    J33 = t327+t115*t324*t371;
    J34 = -t115*t324*t381;
    J35 = rx*t323-(t115*t267*t299*t324)/2.0;
    J40 = -t116*t324*t372;
    J41 = t326+t116*t324*t382;
    J42 = -t323*vy+(t116*t269*t299*t324)/2.0;
    J43 = t116*t324*t371;
    J44 = t327-t116*t324*t381;
    J45 = ry*t323-(t116*t267*t299*t324)/2.0;
    J50 = t393*t405*(t365*(t9+ry*t368+ry*vy)-t362*t366*(t3+t322+rx*t368));
    J51 = t393*t405*(t365*(t3-t4*2.0+t298-(ry*t280*t299)/2.0)+t362*t366*(t9+rx*t370+rx*vx));
    J52 = -t393*t405*(t365*(-t5+t6*2.0+(ry*t269*t299)/2.0)-t362*t366*(t7-t8*2.0+(rx*t269*t299)/2.0));
    J53 = -t393*t405*(t365*(t13+t16+ry*t367)-t362*t366*(t2+rx*t367));
    J54 = t393*t405*(t365*(t2+ry*t369)-t362*t366*(t10+t16+rx*t369));
    J55 = t393*t405*(t365*(rx*rz+(ry*t267*t299)/2.0)+t362*t366*(ry*rz-(rx*t267*t299)/2.0));
  
    return SA[
        J00 J01 J02 J03 J04 J05;
        J10 J11 J12 J13 J14 J15;
        J20 J21 J22 J23 J24 J25;
        J30 J31 J32 J33 J34 J35;
        J40 J41 J42 J43 J44 J45;
        J50 J51 J52 J53 J54 J55
    ]
end