'use strict';
const spawn = require('buffered-spawn');

/**
 * 获取 git log 排除merge, 取最后一条
 *
 * @param   {String} [cwd] git 命令运行时的工作目录
 * @returns {Object} 最后一次提交信息，包含作者名、email、提交时间，提交hash
 */
function runGitLog (cwd) {
	const args = [
		'--no-pager',
		'log',
		'--no-merges',
		'--max-count=1',
		'--format={"hash":"%H","name":"%aN","email":"%aE","time":%at}',
	];
	return spawn('git', args, {
		crossSpawn: false,
		cwd,
	}).then(cp => JSON.parse(cp.stdout.toString()));
}

/**
 * 获取环境变量中的提交信息
 *
 * @returns {Object} 当前提交信息，包含作者名、email、提交时间
 */
function getEvnAuthor () {
	if (/^@?(\d+)(?:\s+[+-]\d+)?$/.test(process.env.GIT_AUTHOR_DATE)) {
		const time = RegExp.$1 - 0;
		const name = process.env.GIT_AUTHOR_NAME;
		const email = process.env.GIT_AUTHOR_EMAIL;
		if (time && name && email) {
			return {
				name: name,
				email: email,
				time: time,
			};
		}
	}
}
const fail = {};
module.exports = function (cwd) {
	return getEvnAuthor() || runGitLog(cwd).catch(error => {
		if (error && error.code === 'ECMDERR') {
			fail[cwd] = error;
		}
	});
};
module.exports.fail = fail;
