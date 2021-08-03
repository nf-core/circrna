'use strict';
/**
 * 过滤错误信息，将明确为其他作者所写的代码、超期的老代码的错误，错误等级标记为警告
 *
 * @param {LintError[]}   errors    已有的错误对象
 * @param {Object}        options   配置
 * @returns {void}
 */
function demoteErrors (errors, options) {
	errors = errors.filter(error => error.blame && !/^0+$/.test(error.blame.hash));
	if (options._expiresTime) {
		errors = errors.filter(error => {
			if (error.blame.author.time < options._expiresTime) {
				demoteError(error);
			} else {
				return error;
			}
		});
	}

	function downMultiline (filter) {
		errors.filter(filter).forEach(demoteError);
	}

	const author = options.author;
	if (author) {
		if (author.email) {
			if (author.email instanceof RegExp) {
				downMultiline(error => !author.email.test(error.blame.author.mail));
			} else {
				downMultiline(error => error.blame.author.mail !== author.email);
			}
		} else if (author.name) {
			if (author.name instanceof RegExp) {
				downMultiline(error => !author.name.test(error.blame.author.name));
			} else {
				downMultiline(error => error.blame.author.name !== author.name);
			}
		}
	}
}

function demoteError (error) {
	error.demote = true;
}

module.exports = demoteErrors;
