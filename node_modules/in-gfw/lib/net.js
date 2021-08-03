"use strict";
const childProcess = require("child_process");
const https = require("https");
const http = require("http");
const url = require("url");
const mem = require("mem");
const syncArgv = "in-gfw-net-sync-argv";

function async (blockedHost, cnHost) {
	let reqs = [];
	function head (config) {
		if (typeof config === "string") {
			if (/^https?:\/\//i.test(config)) {
				config = url.parse(config);
			} else {
				config = {
					hostname: config,
				};
			}
		}
		const request = ((!config.protocol || /^https:$/i.test(config.protocol)) ? https : http).request;
		return new Promise((resolve, reject) => {
			const req = request(Object.assign(config, {
				method: "HEAD",
			}), () => {
				resolve(true);
				abort();
			});

			req.on("error", error => {
				/* istanbul ignore if */
				if (error.errno && /^(?:ETIMEDOUT|ECONNRESET)$/.test(error.errno)) {
					resolve(false);
				} else {
					reject(error);
				}
				abort();
			});
			req.end();
			reqs.push(req);
		});
	}
	function abort () {
		if (reqs) {
			reqs.forEach(req => req.abort());
			reqs = null;
		}
	}
	return Promise.race([
		head(blockedHost || "www.npmjs.com").then(result => !result),
		head(cnHost || "npm.taobao.org"),
	]);
};

function sync (blockedHost, cnHost) {
	const args = [
		__filename,
		syncArgv,
	].concat(
		[
			blockedHost,
			cnHost,
		].filter(Boolean).map(JSON.stringify)
	);
	let result;
	try {
		result = childProcess.execFileSync(
			process.execPath,
			args,
			{
				stdio: "pipe",
				encoding: "utf8",
			}
		);
	} catch (ex) {
		const errorInfo = JSON.parse(ex.stderr);
		const error = new Error(errorInfo.message);
		Object.assign(
			error,
			errorInfo
		);
		throw error;
	}
	return JSON.parse(
		result
	);
}
if ((!process.mainModule || process.mainModule === module) && process.argv[2] === syncArgv) {
	async.apply(this, process.argv.slice(3).map(JSON.parse)).then(result => {
		process.stdout.write(JSON.stringify(result));
		return 0;
	}, error => {
		const errorInfo = {};
		[
			"message",
			"stack",
			"fileName",
			"lineNumber",
			"columnNumber",
		].forEach(prop => {
			errorInfo[prop] = error[prop];
		});
		process.stderr.write(JSON.stringify(Object.assign(errorInfo, error)));
		return 1;
	}).then(code => {
		process.exit(code);
	});
} else {
	module.exports = {
		async: mem(async),
		sync: mem(sync),
	};
}
