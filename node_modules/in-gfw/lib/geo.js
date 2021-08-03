"use strict";
const childProcess = require("child_process");
const regKey = "\\Control Panel\\International\\Geo";
const valueName = "Nation";
const currGeo = "HKCU" + regKey;
const defaultGeo = "HKU\\.DEFAULT" + regKey;

function closeArgsToError (code, stdout, stderr) {
	// code === null when child_process is killed
	if (code) {
		if (stderr) {
			stderr = stderr.replace(/^ERROR:\s*/, "");
		} else {
			stderr = "Exited with status " + code;
		}
		const err = new Error(stderr);
		err.stdout = stdout;
		err.exitStatus = code;
		return err;
	}
	return null;
};
function prepareStream (stream) {
	if (stream == null) {
		return null;
	}
	const buffers = [];
	stream.on("data", (data) => {
		buffers.push(data);
	});
	return buffers;
};

function concatBuffer (buffer) {
	return buffer && Buffer.concat(buffer).toString();
}

function spawn () {
	const args = Array.from(arguments);
	const callback = args.pop();
	const cp = childProcess.spawn.apply(childProcess, args);
	let error;
	cp.on("error", (err) => {
		error = err;
		callback(error);
	});
	let stderr = prepareStream(cp.stderr);
	let stdout = prepareStream(cp.stdout);

	cp.once("exit", code => {
		if (error) {
			return;
		}
		stdout = concatBuffer(stdout);
		stderr = concatBuffer(stderr);
		callback(closeArgsToError(code, stdout, stderr), stdout);
	});
}

function spawnSync () {
	const args = Array.from(arguments);
	const last = args.length - 1;
	const callback = args[last];
	args[last] = {
		encoding: "utf8",
	};
	const cp = childProcess.spawnSync.apply(childProcess, args);
	callback(cp.error || closeArgsToError(cp.status, cp.stdout, cp.stderr), cp.stdout);
}

function geo (spawn) {
	function queryRegFromPowershell (KeyName, valueName, callback) {
		spawn("powershell.exe", [
			"-NoProfile",
			"-ExecutionPolicy",
			"-Command",
			`& {(Get-ItemProperty -Path "Registry::${KeyName}" -Name ${valueName}).${valueName}}`,
		], (error, stdout) => {
			if (stdout && /^(\d+)(?:\r?\n)*$/.test(stdout)) {
				callback(null, RegExp.$1);
			} else {
				callback(error || stdout);
			}
		});
	}
	function queryRegFromRegExe (KeyName, valueName, callback) {
		const args = [
			"QUERY",
			KeyName,
			"/v",
			valueName,
		];
		if (/64$/.test(process.env.PROCESSOR_ARCHITEW6432 || process.arch)) {
			args.push("/reg:64");
		}
		spawn("reg.exe", args, (error, stdout) => {
			if (stdout && /^\s*\S+\s+REG(?:_[A-Z]+)+\s+(.*)$/im.test(stdout)) {
				callback(null, RegExp.$1);
			} else {
				callback(error || stdout);
			}
		});
	}

	function getNation (callback) {
		const opts = [
			queryRegFromRegExe.bind(this, currGeo),
			queryRegFromPowershell.bind(this, currGeo),
			queryRegFromRegExe.bind(this, defaultGeo),
			queryRegFromPowershell.bind(this, defaultGeo),
		];
		function loop () {
			opts.shift()(valueName, (error, nation) => {
				if (error) {
					if (opts.length) {
						loop();
					} else {
						callback(error);
					}
				} else {
					callback(null, nation);
				}
			});
		}
		loop();
	}

	return getNation;
}

module.exports = {
	async: geo(spawn),
	sync: geo(spawnSync),
};
