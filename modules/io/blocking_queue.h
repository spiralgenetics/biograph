
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

// Taken from https://juanchopanzacpp.wordpress.com/2013/02/26/concurrent-queue-c11/ with modifications
 
template <typename T>
class blocking_queue 
{
public:
	// Strong exception safe version of pop
	void pop(T& item)
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		while (queue_.empty())
		{
			cond_.wait(mlock);
		}
		item = queue_.front();
		queue_.pop();
	}

	// Push an item by const ref 
	void push(const T& item)
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		queue_.push(item);
		mlock.unlock();
		cond_.notify_one();
	}

	// Push an item with move semantics
	void push(T&& item) 
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		queue_.push(std::move(item));
		mlock.unlock();
		cond_.notify_one();
	}
private:
	std::queue<T> queue_;
	std::mutex mutex_;
	std::condition_variable cond_;
};


