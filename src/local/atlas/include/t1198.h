const int N = 1198;
static int nloen[N] =
{
  18,
  30,
  32,
  40,
  48,
  60,
  64,
  72,
  80,
  90,
  96,
 100,
 108,
 120,
 120,
 128,
 144,
 144,
 150,
 160,
 160,
 180,
 180,
 180,
 192,
 200,
 200,
 216,
 216,
 240,
 240,
 240,
 240,
 250,
 256,
 270,
 270,
 288,
 288,
 288,
 300,
 300,
 320,
 320,
 320,
 360,
 360,
 360,
 360,
 360,
 360,
 384,
 384,
 384,
 400,
 400,
 400,
 432,
 432,
 432,
 432,
 432,
 450,
 450,
 450,
 480,
 480,
 480,
 480,
 486,
 500,
 500,
 500,
 512,
 512,
 540,
 540,
 540,
 540,
 576,
 576,
 576,
 576,
 576,
 600,
 600,
 600,
 600,
 640,
 640,
 640,
 640,
 640,
 640,
 640,
 648,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 750,
 750,
 750,
 750,
 750,
 768,
 768,
 768,
 800,
 800,
 800,
 800,
 800,
 800,
 810,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 900,
 900,
 900,
 900,
 900,
 900,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 972,
 972,
1000,
1000,
1000,
1000,
1000,
1024,
1024,
1024,
1024,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1280,
1280,
1280,
1280,
1280,
1296,
1296,
1296,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1458,
1458,
1458,
1458,
1500,
1500,
1500,
1500,
1500,
1500,
1500,
1500,
1536,
1536,
1536,
1536,
1536,
1536,
1536,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1620,
1620,
1620,
1620,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1944,
1944,
1944,
1944,
1944,
1944,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2400,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2304,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2250,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2160,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2048,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
2000,
1944,
1944,
1944,
1944,
1944,
1944,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1920,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1800,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1728,
1620,
1620,
1620,
1620,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1600,
1536,
1536,
1536,
1536,
1536,
1536,
1536,
1500,
1500,
1500,
1500,
1500,
1500,
1500,
1500,
1458,
1458,
1458,
1458,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1440,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1350,
1296,
1296,
1296,
1280,
1280,
1280,
1280,
1280,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1250,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1200,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1152,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1080,
1024,
1024,
1024,
1024,
1000,
1000,
1000,
1000,
1000,
 972,
 972,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 960,
 900,
 900,
 900,
 900,
 900,
 900,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 864,
 810,
 800,
 800,
 800,
 800,
 800,
 800,
 768,
 768,
 768,
 750,
 750,
 750,
 750,
 750,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 720,
 648,
 640,
 640,
 640,
 640,
 640,
 640,
 640,
 600,
 600,
 600,
 600,
 576,
 576,
 576,
 576,
 576,
 540,
 540,
 540,
 540,
 512,
 512,
 500,
 500,
 500,
 486,
 480,
 480,
 480,
 480,
 450,
 450,
 450,
 432,
 432,
 432,
 432,
 432,
 400,
 400,
 400,
 384,
 384,
 384,
 360,
 360,
 360,
 360,
 360,
 360,
 320,
 320,
 320,
 300,
 300,
 288,
 288,
 288,
 270,
 270,
 256,
 250,
 240,
 240,
 240,
 240,
 216,
 216,
 200,
 200,
 192,
 180,
 180,
 180,
 160,
 160,
 150,
 144,
 144,
 128,
 120,
 120,
 108,
 100,
  96,
  90,
  80,
  72,
  64,
  60,
  48,
  40,
  32,
  30,
  18
};