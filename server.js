/**
 * MitoRestorer v3.0 — Zero-Dependency Static & Cloud Run Server
 * Supports Google AI Studio Build (*.ai.studio / Cloud Run) and local preview.
 */
const http = require('http');
const fs = require('fs');
const path = require('path');

const PORT = parseInt(process.env.PORT || '8080', 10);
const ROOT_DIR = __dirname;

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
