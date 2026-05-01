import worker from './src/index.js';
const req = new Request('http://localhost/nomogram');
const res = await worker.fetch(req, {}, {});
const html = await res.text();
console.log(html);
