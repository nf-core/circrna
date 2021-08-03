"use strict";
const os = require("./os");
const net = require("./net");

function inGFW (blockedHost, cnHost) {
	return os.async().catch(() => net.async(blockedHost, cnHost));
};
inGFW.sync = function (blockedHost, cnHost) {
	try {
		return os.sync();
	} catch (ex) {
		return net.sync(blockedHost, cnHost);
	}
};
inGFW.os = os.async;
inGFW.osSync = os.sync;
inGFW.net = net.async;
inGFW.netSync = net.sync;

module.exports = inGFW;
