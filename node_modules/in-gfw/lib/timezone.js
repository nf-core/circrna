"use strict";

const fs = require("fs");
const glob = require("glob");
const path = require("path").posix;

function relative (file) {
	return file.replace(/^\/usr\/share\/zoneinfo\//, "");
}

function timezone (fs) {
	function readTimezone (callback) {
		fs.readFile("/etc/timezone", "utf8", (error, data) => {
			if (error) {
				readLocaltime(callback);
			} else {
				callback(null, data.trim());
			}
		});
	}

	function readLocaltime (callback) {
		fs.readlink("/etc/localtime", (error, data) => {
			if (error) {
				findLocaltime(callback);
			} else {
				callback(null, data);
			}
		});
	}

	function lsFiles (pattern, callback) {
		fs.glob(pattern, (error, files) => {
			if (error) {
				callback(error);
			} else {
				callback(null, files);
			}
		});
	}

	function findLocaltime (callback) {
		fs.readFile("/etc/localtime", (error, data) => {
			if (error) {
				callback(error);
			} else {
				lsFiles("/usr/share/zoneinfo/*/*", (error, files) => {
					if (error) {
						callback(error);
					} else {
						findFile(files, data, callback);
					}
				});
			}
		});
	}

	function findFile (files, sourceData, callback) {
		const file = files.shift();
		if (!file) {
			callback(new Error("cannot determine this system's timezone"));
			return;
		}
		fs.stat(file, (error, stats) => {
			if (error || stats.size !== sourceData.length) {
				findFile(files, sourceData, callback);
			} else {
				fs.readFile(file, (error, data) => {
					if (error || !sourceData.equals(data)) {
						findFile(files, sourceData, callback);
					} else {
						callback(null, file);
					}
				});
			}
		});
	}
	return function (callback) {
		readTimezone((error, result) => {
			if (error) {
				callback(error);
			} else {
				result = path.resolve("/usr/share/zoneinfo", result);
				fs.readlink(result, (error, linkString) => {
					if (error) {
						callback(null, relative(result));
					} else {
						callback(null, relative(path.resolve(path.dirname(result), linkString)));
					}
				});
			}
		});
	};
}

const fsSync = {};

function fnSync (fs, fnName) {
	return function () {
		const args = Array.from(arguments);
		const callback = args.pop();
		try {
			callback(null, fs[fnName].apply(fs, args));
		} catch (ex) {
			callback(ex);
		}
	};
}

[
	"readFile",
	"readlink",
	"stat",
].map(fnName => {
	fsSync[fnName] = fnSync(fs, fnName + "Sync");
});

fsSync.glob = fnSync(glob, "sync");

module.exports = {
	sync: timezone(
		fsSync
	),
	async: timezone(
		Object.assign(
			{
				glob,
			},
			fs
		)
	),
};
