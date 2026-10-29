/**
 * MitoRestorer v3.0 — Zero-Dependency Static & Cloud Run Server
 * Supports Google AI Studio Build (*.ai.studio / Cloud Run) and local preview.
 * Automatically proxies /docs/weekly_science_radar_feed.json from GitHub main branch
 * (with 5-minute cache and local disk fallback) so weekly Monday cron updates appear
 * live on *.ai.studio without requiring a manual container republish.
 */
const http = require('http');
const https = require('https');
const fs = require('fs');
const path = require('path');

const PORT = parseInt(process.env.PORT || '8080', 10);
const ROOT_DIR = __dirname;
const GITHUB_FEED_URL = 'https://raw.githubusercontent.com/deesatzed/mitorestorer/main/docs/weekly_science_radar_feed.json';

let cachedRadarFeed = { body: null, ts: 0 };

const MIME_TYPES = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'application/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.cif': 'chemical/x-mmcif; charset=utf-8',
  '.md': 'text/markdown; charset=utf-8',
  '.txt': 'text/plain; charset=utf-8',
  '.ico': 'image/x-icon'
};

function serveLatestRadarFeed(res, localPath) {
  const now = Date.now();
  if (cachedRadarFeed.body && now - cachedRadarFeed.ts < 300000) {
    res.writeHead(200, {
      'Content-Type': 'application/json; charset=utf-8',
      'Cache-Control': 'public, max-age=60',
      'Access-Control-Allow-Origin': '*',
      'X-MitoRestorer-Feed-Source': 'github-main-cache'
    });
    res.end(cachedRadarFeed.body);
    return;
  }

  https.get(GITHUB_FEED_URL, { timeout: 3500 }, (ghRes) => {
    if (ghRes.statusCode !== 200) {
      ghRes.resume();
      serveLocalFallback();
      return;
    }
    let raw = '';
    ghRes.setEncoding('utf8');
    ghRes.on('data', (chunk) => { raw += chunk; });
    ghRes.on('end', () => {
      try {
        JSON.parse(raw);
        cachedRadarFeed = { body: raw, ts: Date.now() };
        res.writeHead(200, {
          'Content-Type': 'application/json; charset=utf-8',
          'Cache-Control': 'public, max-age=60',
          'Access-Control-Allow-Origin': '*',
          'X-MitoRestorer-Feed-Source': 'github-main-live'
        });
        res.end(raw);
      } catch (_) {
        serveLocalFallback();
      }
    });
  }).on('error', () => {
    serveLocalFallback();
  });

  function serveLocalFallback() {
    fs.readFile(localPath, (err, buf) => {
      if (err) {
        res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
        res.end('404 Not Found');
        return;
      }
      res.writeHead(200, {
        'Content-Type': 'application/json; charset=utf-8',
        'Cache-Control': 'public, max-age=60',
        'Access-Control-Allow-Origin': '*',
        'X-MitoRestorer-Feed-Source': 'local-disk-fallback'
      });
      res.end(buf);
    });
  }
}

const server = http.createServer((req, res) => {
  try {
    const parsedUrl = new URL(req.url, `http://${req.headers.host || 'localhost'}`);
    let safePath = decodeURIComponent(parsedUrl.pathname);
    if (safePath === '/' || safePath === '') {
      safePath = '/index.html';
    }

    const fullPath = path.resolve(ROOT_DIR, '.' + safePath);
    if (!fullPath.startsWith(ROOT_DIR)) {
      res.writeHead(403, { 'Content-Type': 'text/plain; charset=utf-8' });
      res.end('403 Forbidden');
      return;
    }

    if (safePath === '/docs/weekly_science_radar_feed.json') {
      serveLatestRadarFeed(res, fullPath);
      return;
    }

    fs.stat(fullPath, (err, stats) => {
      if (err || !stats.isFile()) {
        res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' });
        res.end('404 Not Found: ' + safePath);
        return;
      }

      const ext = path.extname(fullPath).toLowerCase();
      const contentType = MIME_TYPES[ext] || 'application/octet-stream';
      res.writeHead(200, {
        'Content-Type': contentType,
        'Cache-Control': 'public, max-age=300',
        'Access-Control-Allow-Origin': '*'
      });
      fs.createReadStream(fullPath).pipe(res);
    });
  } catch (e) {
    res.writeHead(500, { 'Content-Type': 'text/plain; charset=utf-8' });
    res.end('500 Internal Server Error');
  }
});

server.listen(PORT, '0.0.0.0', () => {
  console.log(`MitoRestorer v3.0 server listening on http://0.0.0.0:${PORT}`);
});
